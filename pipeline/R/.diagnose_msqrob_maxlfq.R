library("dplyr")
library("ggplot2")
library("optparse")
library("readxl")
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "thesis", "genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

RSqM_signif <- read.table(file.path(root_dir, exp_name, "quantification", "RSqM_signif_annotated.tsv"), sep = "\t", header=T)

RSqM_signif_error <- RSqM_signif %>% filter(taxon != "Homo sapiens") %>% arrange(abs(estimate)) %>% filter(abs(estimate) < 0.1)

missed_ecoli <- RSqM_signif_error$Protein.IDs

output_moff <- file.path(root_dir, exp_name, "peptideShaker_out/PSM_reports/output_moff_RAW", "peptide_summary_intensity_moFF_run.tab")
moff <- read.table(output_moff, header = T, sep = "\t")

moff_missed_ecoli <- moff %>% filter(prot %in% missed_ecoli)

peptides <- read.table(file.path(root_dir, exp_name, "paper", "peptides.txt"), sep = "\t", header=T,stringsAsFactors = F)
proteinGroups <- read.table(file.path(root_dir, exp_name, "paper", "proteinGroups.txt"), sep = "\t", header=T, stringsAsFactors = F)
peptides <- peptides %>% select(Sequence, Protein.group.IDs) %>% rename(id = Protein.group.IDs)
peptides$id <- peptides$id %>% as.integer
proteinGroups <- proteinGroups %>% select(id, Protein.IDs)

colnames(proteinGroups)
proteinGroups$id

nrow(moff_missed_ecoli)
peptides_protein <- left_join(peptides, proteinGroups, by = "id") %>% filter(Sequence %in% moff_missed_ecoli$peptide)
source(file.path(root_dir, "scripts", "R", "check_organism.R"))
organism_check <- check_organism(proteins = peptides_protein$Protein.IDs, split=";")
peptides_protein

supplementary <- read_xlsx(file.path(root_dir, exp_name, "mcp.M113.031591-1.xlsx"))

matched <- lapply(moff_missed_ecoli$prot %>% as.character(), function(x) {
  grep(pattern = x, x = supplementary$`Protein IDs`)
  })

matched <- matched %>% lapply(., function(x) ifelse(length(x) == 0, NA, x))

taxonomy <- lapply(matched, function(i) supplementary[[i, "Taxonomy"]])
taxonomy %>% unlist %>% table
