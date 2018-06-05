library("ggplot2")
library("dplyr")
library("tidyr")
library("optparse")
theme_set(theme_bw())

option_list = list(
  make_option(c("--root_dir"), type="character", default="//hest/aoj//thesis/genedata"),
  make_option(c("--exp_name"), type="character", default="thp1"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

input_dir <- ifelse("input_dir" %in% names(opt), opt$input_dir, file.path(root_dir, exp_name, "quantification"))
output_dir <- ifelse("output_dir" %in% names(opt), opt$output, file.path(root_dir, exp_name, "quantification"))


############################################################################
## Read output of custom MaxLFQ
############################################################################

# index <- read.table("protein_ratios.tsv", sep = "\t") %>% .[,ncol(.)]
# Read from Custom LFQ output
protein_intensities <- read.table(file.path(input_dir, "protein_intensities.tsv"), sep = "\t")

############################################################################
## Require intensity in at least 2 of the replicates in each condition
############################################################################
valid_proteins <- (rowSums(protein_intensities[,1:3] != 0) >= 2) & (rowSums(protein_intensities[,4:6] != 0) >= 2)
protein_intensities <- protein_intensities[valid_proteins, ]

############################################################################
## Compute fold change
############################################################################
protein_intensities$log2FC <- log2(rowMeans(protein_intensities[,1:3]) / rowMeans(protein_intensities[,4:6]))
# colnames(protein_intensities) <- c(paste("H", 1:3, sep = ""), paste("L", 1:3, sep = ""), "tax", "fold_change")

############################################################################
## Plot histogram of fold changes segregated by taxonomy
############################################################################

  ## Escherichia coli should have log2FC around log2(3)
  ## Human should remain in log2(1) = 0

  ggplot(data = protein_intensities, aes(x=log2FC)) +
  geom_histogram(alpha=0.5, position="dodge", bins = 50) +
  facet_grid(facets = source ~ .) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-3, 3, 0.5)) +
  theme_gray()
  