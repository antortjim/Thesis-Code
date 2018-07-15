library("MSqRob")
library("dplyr")
library("optparse")
library("MSnbase")

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "thesis", "genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq"),
  make_option(c("--moff_file"), type="character", default="peptide_summary_intensity_moFF_run.tab"),
  make_option(c("--sample_filter"), type="character", default=""),
  make_option(c("--experiment_contrasts"), type="character", default="conditionH-conditionL",
              help = "String with structure:
              condition condition1 - condition condition2
              with no spaces (spaces are shown for illustrative purposes"),
  make_option(c("--save_model"), action="store_true", default=TRUE, help="Do you want to save the model?"),
  make_option(c("--suffix"), type="character", help="Suffix to append to the RSqM_signif file", default=""),
  make_option(c("--fraction_normalized"), action="store_true", default=TRUE,
              help="Are you passing a file that has been normalized using MaxLFQ Levenberg-Marquandt minimisation?"),
  make_option(c("--normalisation"), default="quantiles",
              help="Should normalisation be performed? Pass \"none\" for no quantilsation"),
 
  make_option(c("--smallest_unique_groups"), action="store_true", default=TRUE,
              help="Should protein groups be simplified?")
 
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
moff_file <- opt$moff_file
sample_filter <- opt$sample_filter
experiment_contrasts <- opt$experiment_contrasts
experiment_contrasts <- ifelse(exp_name == "thp1", "conditionpLPS_pellet-conditionmLPS_pellet", experiment_contrasts)
save_model <- opt$save_model
suffix <- opt$suffix
fraction_normalized <- opt$fraction_normalized
smallest_unique_groups <- opt$smallest_unique_groups
normalisation <- opt$normalisation
moff_file <- ifelse(fraction_normalized, "peptide_summary_intensity_moFF_run_fraction_normalized.tab", "peptide_summary_intensity_moFF_run.tab")
fraction_normalized <- ifelse(exp_name == "thp1", F, fraction_normalized)

output_moff <- file.path(root_dir, exp_name, "peptideShaker_out/PSM_reports/output_moff_RAW")

print("Importing MSnSet")
peptides_msnset <- import2MSnSet(file = file.path(output_moff, moff_file),
                          filetype = "moFF")

exp_annot <- read.table(file.path(root_dir, exp_name, "data", "experimental_design.tsv"), header = T) %>%
  select(-Group)

colnames(exp_annot) <- c("run", "fraction", "condition", "rep")
if(!sample_filter == "") {
  exp_annot <- exp_annot[grep(pattern = sample_filter, exp_annot$Experiment),]
}

if(fraction_normalized) {
  exp_annot <- exp_annot %>% filter(fraction==1)
  exp_annot$run <- factor(paste0(exp_annot$condition, exp_annot$rep))
}

print("Preprocessing MSnSet")
# Takes the log2 of the intensities, performs quantile.normalisation,
# adds the experimental annotation, removes peptides identified in only 1 sample
# solves the protein inference problem by taking the smallestUniqueGroups i.e
# if a protein appears in several protein groups, only the smallest one is kept
peptides_msnset_processed <- preprocess_MSnSet(MSnSet = peptides_msnset,
                                            accession = "prot",
                                            exp_annotation = exp_annot,
                                            logtransform = TRUE, base = 2,
                                            normalisation = normalisation,
                                            smallestUniqueGroups = smallest_unique_groups, split = ", ",
                                            useful_properties = c("prot", "peptide"),
                                            minIdentified = 2)

print("Compiling protdata object")
proteins_protdata <- MSnSet2protdata(peptides_msnset_processed, accession="prot")

#Fixed effects
fixed <- c("condition")

#Random effects, for label-free data, it is best to always keep "Sequence"
random <- c("peptide")
if(fraction_normalized) {random <- c(random, "fraction"); print("Adding fraction as random effect")}



#To which folder do you want to export the Excel file(s) and the saved model file (if you chose to do so)? Please do not use a trailing "/" here!
export_folder <- file.path(root_dir, exp_name, "quantification")

#Construct the contrast matrix L for testing on the fold changes of interest (i.e. our research hypotheses)
L <- makeContrast(contrasts=experiment_contrasts,
                  levels=strsplit(experiment_contrasts, split = "-") %>% unlist)

#Set the significance threshold (default: 5% FDR)
# FDRlevel <- 0.05

print("Fitting ridge regression model with Huber weights and empirical Bayes estimation of variance")

system.time(model_RR <- fit.model(protdata=proteins_protdata,
                                  response="quant_value",
                                  fixed=fixed, random=random,
                                  shrinkage.fixed=NULL,
                                  weights="Huber",
                                  squeezeVar=TRUE))

print(paste0("Model fit for n proteins: ", length(model_RR)))


print("Performing hypothesis testing")
RSqM <- test.protLMcontrast(model_RR, L)
print("Multiple testing correction")
RSqM_adjust <- prot.p.adjust(RSqM)
RSqM_signif <- prot.signif(RSqM_adjust)
RSqM_signif$Protein.IDs <- rownames(RSqM_signif)
print("DONE")
print(paste0("Saving results to ", export_folder))


write.table(x = RSqM_adjust, file = file.path(export_folder, paste0("RSqM_adjust", suffix, ".tsv")),
            sep = "\t", col.names = T, row.names = F, quote = F)


write.table(x = RSqM_signif, file = file.path(export_folder, paste0("RSqM_signif", suffix, ".tsv")),
            sep = "\t", col.names = T, row.names = F, quote = F)


#If you chose to save the model, save it
if(isTRUE(save_model)){
  print("Saving model")
  result_files <- list()
  result_files$proteins <- proteins_protdata
  result_files$models <- model_RR
  result_files$levelOptions <- rownames(L)
  result_files$fixed <- fixed
  result_files$random <- random
  saves_MSqRob(result_files, file=file.path(export_folder,"model.RDatas"), overwrite=TRUE)
  print("Model saved")
}
