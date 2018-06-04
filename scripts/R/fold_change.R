library("ggplot2")
library("dplyr")
library("tidyr")
# theme_set(theme_bw())
# source("load_maxlfq_data.R")

############################################################################
## Read output of custom MaxLFQ
############################################################################

# index <- read.table("protein_ratios.tsv", sep = "\t") %>% .[,ncol(.)]
# Read from Custom LFQ output
protein_intensities <- read.table("protein_intensities.tsv", sep = "\t")
log2ratioLFQ <- read.table("log2ratioLFQ.tsv", sep = "\t")[,1]
taxonomy <- read.table("Taxonomy.tsv", sep = "\t", stringsAsFactors = F)[,1]

# browser()
protein_intensities$taxonomy <- taxonomy
# protein_intensities$index <- index

############################################################################
## Require intensity in at least 2 of the replicates in each condition
############################################################################
valid_proteins <- (rowSums(protein_intensities[,1:3] != 0) >= 2) & (rowSums(protein_intensities[,4:6] != 0) >= 2)
protein_intensities <- protein_intensities[valid_proteins, ]

############################################################################
## Compute fold change
############################################################################
protein_intensities$custom <- log2(rowMeans(protein_intensities[,1:3]) / rowMeans(protein_intensities[,4:6]))
# colnames(protein_intensities) <- c(paste("H", 1:3, sep = ""), paste("L", 1:3, sep = ""), "tax", "fold_change")
protein_intensities$MaxLFQ <- log2ratioLFQ[valid_proteins]

protein_intensities <- gather(protein_intensities, source, log2FC, -(V1:taxonomy))

protein_intensities %>% group_by(source, taxonomy) %>% summarise(count = n())
protein_intensities %>% group_by(source, taxonomy) %>% summarise(meanlog2FC = mean(log2FC))
protein_intensities %>% group_by(source, taxonomy) %>% summarise(stdlog2FC = sd(log2FC))

############################################################################
## Plot histogram of fold changes segregated by taxonomy
############################################################################

  ## Escherichia coli should have log2FC around log2(3)
  ## Human should remain in log2(1) = 0

  ggplot(data = protein_intensities, aes(x=log2FC, fill=taxonomy)) +
  geom_histogram(alpha=0.5, position="dodge", bins = 50) +
  facet_grid(facets = source ~ .) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-3, 3, 0.5)) +
  theme_gray()

# Average log2FC
protein_intensities %>% group_by(source, taxonomy) %>% summarise(average_log2FC = mean(log2FC))

# Counts of protein groups per taxonomy
protein_intensities %>% group_by(source, taxonomy) %>% summarise(count = n())
