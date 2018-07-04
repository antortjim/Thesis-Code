library("ggplot2")
library("dplyr")
library("tidyr")
library("optparse")
library("readxl")
theme_set(theme_bw())
source(file.path(root_dir, "scripts", "R", "map_to_organism.R"))

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character"),
  make_option(c("--plots_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

input_dir <-  file.path(root_dir, exp_name, ifelse("input_dir" %in% names(opt), opt$input_dir, "quantification"))
output_dir <- file.path(root_dir, exp_name, ifelse("output_dir" %in% names(opt), opt$output_dir, "quantification"))
plots_dir <- file.path(root_dir, exp_name, ifelse("output_dir" %in% names(opt), opt$plots_dir, "plots"))


############################################################################
## Read protein quantification data
############################################################################

###
msnbase_q <- read.table(file = file.path(root_dir, exp_name, "quantification", "msnbase_quantification.tsv"), sep = "\t", header=T)
###
lfq_q <- read.table(file = file.path(root_dir, exp_name, "quantification", "protein_intensities_own_pipeline.tsv"), sep = "\t", header=T)
rownames(lfq_q) <- lfq_q[,1]
lfq_q <- lfq_q[,-1]
###
hybrid <- read.table(file = file.path(root_dir, exp_name, "quantification", "protein_intensities_maxlfq_paper.tsv"), sep = "\t", header=T)
rownames(hybrid) <- hybrid[,1]
hybrid <- hybrid[,-1]
###
extract_taxonomy <- function(quantification, split = ", ") {
  Protein.IDs <- rownames(quantification)
  taxonomy <- lapply(Protein.IDs, function(x) check_id_taxonomy(x, id_mapping, split)) %>% unlist
  inferred_taxonomy <- data.frame(Protein.IDs = Protein.IDs, Taxonomy = taxonomy)
  return(inferred_taxonomy)
}

############################################################################
## Require intensity in at least 2 of the replicates in each condition
############################################################################
process_intensities <- function(quantification, inferred_taxonomy) {
  Protein.IDs <- rownames(quantification)
  colnames(quantification) <- gsub("_", "", colnames(quantification))
  
  protein_intensities <- quantification
  valid_proteins <- (rowSums(protein_intensities[,1:3] != 0) >= 2) & (rowSums(protein_intensities[,4:6] != 0) >= 2)
  protein_intensities <- protein_intensities[valid_proteins, ]
  Protein.IDs <- Protein.IDs[valid_proteins]
  
  ############################################################################
  ## Compute fold change
  ############################################################################
  
  protein_intensities$log2FC <- log2(rowMeans(protein_intensities[,1:3]) / rowMeans(protein_intensities[,4:6]))
  protein_intensities$Protein.IDs <- Protein.IDs
  
  ## Join taxonomy
  protein_intensities <- left_join(protein_intensities, inferred_taxonomy, by = "Protein.IDs")
}

inferred_taxonomy <- extract_taxonomy(lfq_q)
lfq_q <- process_intensities(lfq_q, inferred_taxonomy)
inferred_taxonomy <- extract_taxonomy(msnbase_q, split=", ")
msnbase_q <- process_intensities(msnbase_q, inferred_taxonomy)


###############################
## Read the published output
###############################

maxlfq_output <- read.table(file.path(root_dir, "maxlfq", "paper", "proteinGroups.txt"), sep = "\t", header=T)
maxlfq_columns <- "LFQ.intensity."
maxlfq_output_filtered <- maxlfq_output[maxlfq_output$Protein.IDs %in% Protein.IDs, c("Protein.IDs", grep(pattern = maxlfq_columns, colnames(maxlfq_output), value = T))]
maxlfq_output_filtered$log2FC <- log2(rowMeans(maxlfq_output_filtered[,2:4]) / rowMeans(maxlfq_output_filtered[,5:7]))
colnames(maxlfq_output_filtered) <- colnames(maxlfq_output_filtered) %>% gsub(pattern = maxlfq_columns, replacement = "")
published_taxonomy <- read.table(file.path(input_dir, "Taxonomy.tsv"), sep = "\t", header = T, stringsAsFactors = F)
maxlfq_output_filtered <- left_join(maxlfq_output_filtered, published_taxonomy, by = "Protein.IDs")


hybrid_intensities <- process_intensities(hybrid, published_taxonomy)


plot_data <- rbind(
  cbind(maxlfq_output_filtered %>% select(Protein.IDs, log2FC, H1:L3, Taxonomy), source = "maxlfq"),
  # cbind(lfq_q %>% select(Protein.IDs, log2FC, H1:L3, Taxonomy), source = "pipeline_lfq"),
  cbind(msnbase_q %>% select(Protein.IDs, log2FC, H1:L3, Taxonomy), source = "pipeline_msnbase"),
  cbind(hybrid_intensities %>% select(Protein.IDs, log2FC, H1:L3, Taxonomy), source = "hybrid")
)

plot_data %>% group_by(source) %>% summarise(sum(is.na(Taxonomy)))

# colnames(protein_intensities) <- c(paste("H", 1:3, sep = ""), paste("L", 1:3, sep = ""), "tax", "fold_change")

############################################################################
## Plot histogram of fold changes segregated by taxonomy
############################################################################

  ## Escherichia coli should have log2FC around log2(3)
  ## Human should remain in log2(1) = 0

  ggplot(data = plot_data %>% filter(!is.na(Taxonomy)), aes(x=log2FC)) +
  geom_histogram(alpha=0.5, position="dodge", bins = 50, aes(fill = Taxonomy)) +
  facet_grid(facets = source ~ .) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-3, 3, 0.5)) +
  theme_gray()

ggsave(filename = file.path(plots_dir, "histogram_maxlfq_custom.png")) 
