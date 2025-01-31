library("MSqRob")
library("dplyr")
library("optparse")
library("MSnbase")
library("latex2exp")


home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj",
                   ifelse(Sys.info()["user"] == "aoj", "/z/home/aoj",
                          "/home/antortjim/"))


data_dir <- ifelse(Sys.info()["user"] == "aoj", "thesis/genedata/maxlfq/",
                   "MEGA/Master/Thesis/Code/scripts/data")

exp_dir <-  ifelse(Sys.info()["sysname"] == "Windows",
                   file.path(root_dir, exp_name),
                   file.path(home_dir, data_dir))

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
fraction_normalized <- ifelse(exp_name == "maxlfq", FALSE, opt$fraction_normalized)
smallest_unique_groups <- opt$smallest_unique_groups
normalisation <- ifelse(opt$normalisation != "maxlfq", opt$normalisation, "none")
normalisation <- opt$normalisation

moff_file <- ifelse(fraction_normalized, "peptide_summary_intensity_moFF_run_fraction_normalized.tab", "peptide_summary_intensity_moFF_run.tab")
fraction_normalized <- ifelse(exp_name == "thp1", F, fraction_normalized)
#To which folder do you want to export the Excel file(s) and the saved model file (if you chose to do so)? Please do not use a trailing "/" here!
export_folder <-  ifelse(Sys.info()["sysname"] == "Windows",
                         file.path(exp_dir, "quantification"),
                         exp_dir
)


output_moff <- ifelse(Sys.info()["sysname"] == "Windows",
                   file.path(exp_dir, "peptideShaker_out/PSM_reports/output_moff_RAW"),
                   file.path(exp_dir, "output_moFF_RAW")
                   )

exp_annot <- read.table(file.path(
  ifelse(Sys.info()["sysname"] == "Windows",
         file.path(root_dir, exp_name, "data"),
         exp_dir),
  "experimental_design.tsv"), header = T) %>%
  select(-Group)

colnames(exp_annot) <- c("run", "fraction", "condition", "rep")
if(!sample_filter == "") {
  exp_annot <- exp_annot[grep(pattern = sample_filter, exp_annot$Experiment),]
}

if(fraction_normalized) {
  exp_annot <- exp_annot %>% filter(fraction==1)
  exp_annot$run <- factor(paste0(exp_annot$condition, exp_annot$rep))
}


print("Importing MSnSet")
  
#### EXPERIMENT
moff_data <- read.table(file.path(output_moff, moff_file), header=T, sep = "\t")

apex_data <- moff_data[,-(1:2)]
colnames(apex_data) <- colnames(apex_data) %>% gsub(pattern = "sumIntensity_", replacement = "")

conds <- c("H", "L")
reps <- exp_annot$rep %>% unique %>% as.integer

max_aggregation <- data.frame(peptide = moff_data[,1], prot = moff_data[,2])
for (cond in conds) {
  for (repl in reps) {
    sample_name <- paste0(cond, repl)
    print(sample_name)
    max_aggregation[,sample_name] <- apply(apex_data[,(exp_annot %>% filter(condition==cond, rep == repl) %>% .$run)], 1, max)
  }
}
colnames(max_aggregation)[-(1:2)] <- paste0("sumIntensity_", colnames(max_aggregation)[-(1:2)])
moff_file <- "peptide_summary_intensity_moFF_run_max.tab"
suffix <- "_max"

exp_annot <- exp_annot %>% filter(fraction == 1) %>% select(condition, rep)
exp_annot$run <- factor(paste0(exp_annot$condition, exp_annot$rep))

write.table(x = max_aggregation, file.path(output_moff, moff_file), col.names = T, sep = "\t", quote=F, row.names=F)
# write.table(x = max_aggregation, file.path(moff_file), col.names = T, sep = "\t", quote=F, row.names=F)

#### EXPERIMENT
peptides_msnset <- import2MSnSet(file = file.path(output_moff, moff_file),
                          filetype = "moFF")


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


MSqRob_quantification <- function(peptides_msnset_processed, fraction_normalized=FALSE, experiment_contrasts, save_model=FALSE) {
  print("Compiling protdata object")
  proteins_protdata <- MSnSet2protdata(peptides_msnset_processed, accession="prot") & system('notify-send "Done"')
  proteins_protdata_none <- MSnSet2protdata(peptides_msnset_processed_none, accession="prot") & system('notify-send "Done"')
  
  #Fixed effects
  fixed <- c("condition")
  
  #Random effects, for label-free data, it is best to always keep "Sequence"
  random <- c("peptide")
  if(!fraction_normalized) {random <- c(random, "fraction"); print("Adding fraction as random effect")}
  
  
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
  RSqM <- test.protLMcontrast(model_RR, L) & system('notify-send "Done"')
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
  return(RSqM_signif)
  
}

RSqM_signif <- MSqRob_quantification(peptides_msnset_processed = peptides_msnset_processed, experiment_contrasts = experiment_contrasts)




