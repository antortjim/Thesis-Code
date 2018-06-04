library("dplyr")
library("optparse")

option_list = list(
  make_option(c("--root_dir"), type="character", default="//hest/aoj//thesis/genedata"),
  make_option(c("--exp_name"), type="character", default="thp1_test"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

input_dir <- ifelse("input_dir" %in% names(opt), opt$input_dir, file.path(root_dir, exp_name, "quantification"))
output_dir <- ifelse("output_dir" %in% names(opt), opt$output, file.path(root_dir, exp_name, "quantification"))

######################################################################
## Read data
######################################################################

intensities <- read.table(file = file.path(input_dir, "intensities.tsv"), sep = "\t", header=T, check.names = F)
combinations <- read.table(file = file.path(input_dir, "combinations.tsv"), sep = "\t") %>% as.matrix %>% .[-1,]
combinations_names <- read.table(file = file.path(input_dir, "combinations_names.tsv"), sep = "\t", stringsAsFactors = F)[,1]
peptides <- read.table(file = file.path(input_dir, "peptides.tsv"), sep = "\t", header=T, stringsAsFactors = F) %>% select(Protein.IDs, Sequence)
Protein.IDs <- read.table(file = file.path(input_dir, "Protein.IDs.tsv"), sep = "\t", header=T)[,1]

######################################################################
## Compute peptide ratios from peptide intensity
######################################################################

peptide_ratios <- combinations %>% apply(X = ., MARGIN = 2, FUN = function(x) intensities[, x[1]] / intensities[, x[2]])
colnames(peptide_ratios) <- combinations_names
rownames(peptide_ratios) <- peptides$Sequence
# if a peptide ratio is infinite, the denominator was 0
# if a peptide ratio is NaN, the denominator and the numerator were 0
# if a peptide ratio is 0, the numerator was 0
peptide_ratios[is.na(peptide_ratios) | is.infinite(peptide_ratios)] <- 0

######################################################################
## Compute protein ratios from peptide ratios
######################################################################

protein_ratios <- matrix(NA, nrow = length(Protein.IDs), ncol = length(combinations_names))
colnames(protein_ratios) <-  combinations_names
rownames(protein_ratios) <- Protein.IDs

minimal_shared_peptides <- 2
i <- 1
for (protein in Protein.IDs) {
  # Get the peptides associated with protein protein
  protein_pos <- which(peptides$Protein.IDs == protein)
  # If the number of associated peptides is NOT below the threshold, 
  if(!length(protein_pos) < minimal_shared_peptides) {
    
    # Fetch the peptide ratios for the current protein
    protein_peptide_ratios <- peptide_ratios[protein_pos,] %>% as.matrix()
    
    # The protein ratio for comparison k is the median of the peptide ratios in comparison k 
    r <-  protein_peptide_ratios %>% apply(X = ., MARGIN = 2, FUN = function(x) median(x, na.rm = T))
    # How many peptides support each median? This number goes from 0 to the total number of peptides for this protein
    # It's the count of non-zero peptides
    number_peptides_supporting_median <- length(protein_pos) - (protein_peptide_ratios %>% apply(X = ., MARGIN = 2, FUN = function(x) x==0) %>% colSums)
    r[number_peptides_supporting_median < minimal_shared_peptides] <- 0
  # else, set the protein ratios equal to 0
  } else {
    r <- rep(0, length(combinations_names))
  }
  
  # Set the ratios of the ith protein in the ith row of the matrix
  protein_ratios[i,] <- r
  
  i <- i + 1
}

######################################################################
## Visualize distribution of the number of ratios equal to 0
## (max equal to the number of combinations, i.e 15) 
######################################################################
x <- rowSums(protein_ratios == 0)
table(x)
hist(x)
# Hopefully most proteins have 0 ratios = 0, i.e the histogram should have a high bin at 0 and very low bins elsewhere
# and lower the more to the right we move in the x axis

######################################################################
## Export the protein ratios to a file
######################################################################
write.table(x = protein_ratios, file = file.path(output_dir, "protein_ratios.tsv"), row.names = T, col.names = T, sep = "\t")
