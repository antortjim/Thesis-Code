library("dplyr")
library("ggplot2")
library("gProfileR")
library("optparse")

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "thesis", "genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
source(file.path(root_dir, "scripts", "R", "check_organism.R"))

RSqM_signif_file <- file.path(root_dir, exp_name, "quantification", "RSqM_signif.tsv")
RSqM_signif <- read.table(file = RSqM_signif_file, header = T, sep = "\t")

organism_check <- check_organism(as.character(RSqM_signif$Protein.IDs), split = ", ")
RSqM_signif$taxon <- organism_check[[2]]

plots_dir <- file.path(root_dir, exp_name, "plots")

volcano_plot <- ggplot(data = RSqM_signif %>% filter(!is.na(taxon)), aes(x = estimate, y = -log10(qval))) +
  # geom_point() +
  geom_point(aes(col = taxon)) +
  geom_vline(xintercept = c(-1, 1)) +
  geom_hline(yintercept = -log10(0.05))

ggsave(filename = file.path(plots_dir, "volcano_plot_moff.png"), plot = volcano_plot)

ggplot(data = RSqM_signif %>% filter(!is.na(taxon)), aes(x = estimate, fill=taxon)) +
  geom_histogram(bins = 60, position="dodge")  +
  scale_x_continuous(breaks = -6:6)

up_prots <- RSqM_signif %>% filter(signif) %>% filter(estimate > 1) %>% .$Protein.IDs %>% strsplit(split = ", ") %>% unlist
down_prots <- RSqM_signif %>% filter(signif) %>% filter(estimate < -1) %>% .$Protein.IDs %>% strsplit(split = ", ") %>% unlist
stable_prots <- RSqM_signif %>% filter(signif) %>% filter(abs(estimate) < 1) %>% .$Protein.IDs %>% strsplit(split = ", ") %>% unlist
dap_prots <- c(up_prots, down_prots)
na_prots <- RSqM_signif %>% filter(is.na(estimate)) %>% .$Protein.IDs

proteinsTHP1[na_prots[1]]
lapply(1:length(na_prots), function(i) proteinsTHP1[na_prots[i]]@data[[1]]$Experiment %>%
         as.character %>% unique)

RSqM_signif %>% filter(estimate == min(estimate, na.rm = T))

gsea <- gprofiler(query = c(dap_prots, stable_prots), organism = "hsapiens")
gsea <- gprofiler(query = c(na_prots), organism = "hsapiens")