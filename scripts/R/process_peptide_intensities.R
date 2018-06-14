library("dplyr")
library("optparse")
library("ggplot2")
library("preprocessCore")

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character"),
  make_option(c("--pipeline"), type="logical", default = F)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
pipeline <- opt$pipeline
# if pipeline false, the information in peptides.txt is read
pipeline <- F
# if pipeline true, the information from the xic normalization step in the pipeline is read
pipeline <- T

input_dir <- file.path(root_dir, exp_name, ifelse("input_dir" %in% names(opt), opt$input_dir, "quantification"))
output_dir <- file.path(root_dir, exp_name, ifelse("output_dir" %in% names(opt), opt$output, "quantification"))

######################################################################
## Read data
######################################################################
experimental_design <- read.table(file = file.path(root_dir, exp_name, "data", "experimental_design.tsv"), header=T)
fractions_available <- F
if(any(experimental_design$Fraction != 1)) fractions_available <- T


## A from MaxQuant paper published results
###############################################################

if(!pipeline) {
  # The intensities file from MaxLFQ paper was generated using the script: /thesis/genedata/maxlfq/paper/load_maxlfq_data.R
  # This loads the protein intensities prior to LFQ processing (least squares)
  intensities <- read.table(file = file.path(input_dir, "maxlfq_peptide_intensities.tsv"), sep = "\t", header=T, check.names = F)
  # The protein ids associated to each peptide i.e the protein column in the PSM report
  peptides_Protein.IDs <- read.table(file = file.path(input_dir, "peptides_Protein.IDs.txt"), check.names = F, stringsAsFactors = F)[,1]
  # The sequence associated to each peptide i.e the sequence column in the PSM report
  peptides_sequence <- read.table(file = file.path(input_dir, "peptides_sequence.txt"), check.names = F, stringsAsFactors = F)[,1]
  
  combinations <- read.table(file = file.path(input_dir, "combinations.tsv"), sep = "\t") %>% as.matrix %>% .[-1,]
  combinations_names <- read.table(file = file.path(input_dir, "combinations_names.tsv"), sep = "\t", stringsAsFactors = F)[,1]
  
  suffix <- "_maxlfq_paper"
  # Protein.IDs <- read.table(file = file.path(input_dir, "Protein.IDs.tsv"), sep = "\t", header=F, stringsAsFactors = F)[,1]
  Protein.IDs <- unique(peptides_Protein.IDs)
} else {
  # read from xic_norm
  intensities <- read.table(file = file.path(root_dir, exp_name, "peptideShaker_out", "PSM_reports",
                                             "output_moFF_raw", "peptide_summary_intensity_moFF_run_xic_norm.tab"),
                            sep = "\t", header=T)
  colnames(intensities) <- colnames(intensities) %>% gsub("_", "", .)
  peptides_Protein.IDs <- intensities$Protein.IDs
  peptides_sequence <- intensities$Sequence
  combinations <- combn(x = colnames(intensities)[-c(1,2)], m = 2)
  combinations_names <- combinations %>% apply(., 2, function(x) paste(x, collapse = "/"))
  suffix <- "_own_pipeline"
  
}

######################################################################
## Preprocess intensities: take the log2 and quantile.normalize()
######################################################################

conditions <- paste0(experimental_design$Experiment, experimental_design$Replicate) %>% unique %>% sort

intensities_matrix <- as.matrix(intensities[conditions])
if(pipeline) {
intensities_matrix <- log2(intensities_matrix)
intensities_matrix[is.infinite(intensities_matrix)] <- NA
}

log2_qnorm_intensities <- normalize.quantiles(intensities_matrix)
# log2_intensities[is.na(log2_intensities)] <- 0

intensities[conditions] <- log2_qnorm_intensities

# intensities[conditions] <- log2_intensities
head(log2_qnorm_intensities)
log2_qnorm_intensities <- log2_qnorm_intensities %>% as.data.frame
colnames(log2_qnorm_intensities) <- conditions

ggplot(data = log2_qnorm_intensities %>% gather(),
       aes(y = value, fill = key, x = key)) + geom_boxplot()
ggsave(file.path(root_dir, exp_name, "plots", "intensity_normalization_across_samples.png"))

write.table(x = intensities, file = file.path(root_dir, exp_name, "peptideShaker_out", "PSM_reports", "output_moFF_raw",
                                              paste0("peptide_summary_intensity_qnorm", suffix, ".tsv")), row.names = T, col.names = T, sep = "\t", quote = F)


# ggplot(data = select(intensities, H1:L3) %>% gather() %>% filter(value != 0),
#        aes(x = value, col = key)) + geom_density()

######################################################################
## Compute peptide ratios from peptide intensity
######################################################################

peptide_ratios <- combinations %>% apply(X = ., MARGIN = 2, FUN = function(x) intensities[, x[1]] / intensities[, x[2]])
colnames(peptide_ratios) <- combinations_names
rownames(peptide_ratios) <- peptides_sequence
# if a peptide ratio is infinite, the denominator was 0
# if a peptide ratio is NaN, the denominator and the numerator were 0
# if a peptide ratio is 0, the numerator was 0
peptide_ratios[is.na(peptide_ratios) | is.infinite(peptide_ratios)] <- 0


######################################################################
## Peptides per protein
######################################################################
# peptides %>% group_by(Protein.IDs) %>% summarise(count = n()) %>% ggplot(aes(x = count)) + geom_histogram(bins = 40)
# peptides %>% group_by(Protein.IDs) %>% summarise(count = n()) %>% .$count %>% table %>% .[-1] %>% sum


######################################################################
## Compute protein ratios from peptide ratios
######################################################################

protein_ratios <- matrix(NA, nrow = length(Protein.IDs), ncol = length(combinations_names))
colnames(protein_ratios) <-  combinations_names
rownames(protein_ratios) <- Protein.IDs
status <- character(length = length(Protein.IDs))
supporting <- integer(length = length(Protein.IDs))
local_supporting <- matrix(nrow = length(Protein.IDs), ncol = ncol(combinations))
minimal_shared_peptides <- 2
i <- 1
for (protein in Protein.IDs) {
  print(i)
  # Get the peptides associated with protein protein
  protein_pos <- which(peptides_Protein.IDs == protein)
  supporting[i] <- length(protein_pos)
  # If the number of associated peptides is NOT below the threshold, 
  if(!supporting[i] < minimal_shared_peptides) {
    
    # Fetch the peptide ratios for the current protein
    protein_peptide_ratios <- peptide_ratios[protein_pos,] %>% as.matrix()
    # if("GYPHWPAR" %in% rownames(protein_peptide_ratios)) browser()
    
    # The protein ratio for comparison k is the median of the peptide ratios in comparison k 
    r <-  protein_peptide_ratios %>% apply(X = ., MARGIN = 2, FUN = function(x) median(x, na.rm = T))
    # How many peptides support each median? This number goes from 0 to the total number of peptides for this protein
    # It's the count of non-zero peptides
    number_peptides_supporting_median <- length(protein_pos) - (protein_peptide_ratios %>% apply(X = ., MARGIN = 2, FUN = function(x) x==0) %>% colSums)
    local_supporting[i,] <- number_peptides_supporting_median
    r[number_peptides_supporting_median < minimal_shared_peptides] <- 0
    status[i] <- "SUPPORT"
    # else, set the protein ratios equal to 0
  } else {
    r <- rep(0, length(combinations_names))
    status[i] <- "FRAGILE"
  }
  
  # Set the ratios of the ith protein in the ith row of the matrix
  protein_ratios[i,] <- r
  
  i <- i + 1
}

supporting %>% table


######################################################################
## Visualize distribution of the number of ratios equal to 0
## (max equal to the number of combinations, i.e 15) 
######################################################################
x <- rowSums(protein_ratios == 0)
table(x == ncol(protein_ratios))
table(x)
hist(x, breaks = 30)
# Hopefully most proteins have 0 ratios = 0, i.e the histogram should have a high bin at 0 and very low bins elsewhere
# and lower the more to the right we move in the x axis
# If most proteins have 0 ratios = 0 means that the aggregation to proteins didnt go well

######################################################################
## Export the protein ratios to a file
######################################################################
write.table(x = protein_ratios, file = file.path(output_dir, paste0("protein_ratios", suffix, ".tsv")), row.names = T, col.names = T, sep = "\t", quote = F)
# Call LFQ.py with this file
# or msnbase
write.table(x = rownames(protein_ratios), file = file.path(output_dir, paste0("protein_ratios.IDs", suffix, ".txt")), row.names = F, col.names = F, sep = "\t", quote = F)
