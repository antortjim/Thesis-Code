library(ggplot2); library(dplyr); library(tidyr); library(ggalt); library(cowplot)
source("plots.R")
bayesquant <- read.table("data/bayesquant_res.tsv", header=T, sep = "\t")
colnames(bayesquant)[1] <- "protein"
print(colnames(bayesquant))
plots_dir <- "../../Report/plots/"

p1 <- make_dumbbell_plot(bayesquant[bayesquant$n_peptides %in% c(2,3,4,6,7,10),], n_pep = NULL, n=5,ncol=2,test=T) 
p1
ggsave(file.path(plots_dir, "performance.eps"), height = 7, width=5, plot = p1)
ggsave(file.path(plots_dir, "performance.png"), height = 7, width=5, plot = p1)
p1 <- make_dumbbell_plot(bayesquant[bayesquant$n_peptides %in% c(2,3,4,6,7,10),], n_pep = NULL, n=5,ncol=3,test=T) 
p1
ggsave(file.path(plots_dir, "performance_horizontal.eps"), height = 5, width=7, plot = p1)
ggsave(file.path(plots_dir, "performance_horizontal.png"), height = 5, width=7, plot = p1)

