library("dplyr")
library("ggplot2")
library("gProfileR")
library("optparse")
source("read_quant.R")
source("plots")
library(latex2exp)
theme_set(theme_bw(base_size=15))
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "thesis", "genedata")),
  make_option(c("--exp_name"), type="character", default="thp1"),
  make_option(c("--output_dir"), type="character", default="."),
  make_option(c("--suffix"), type="character", default=""),
  make_option(c("--quant_file"), type="character", default=""),
  make_option(c("--rope"), default=.4)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(interactive()) opt <- list(rope = 0.4, quant_file = "../thp1/bayesquant_res.tsv",
                              output_dir = "../thp1", suffix="_bq")

root_dir <- opt$root_dir
exp_name <- opt$exp_name
suffix <- opt$suffix
plots_dir <- file.path(opt$output_dir, "plots")
rope <- opt$rope


quant <- read_quant(filename = opt$quant_file, filetype = "BayesQuant")
nrow(quant)

out_rope_log <- (quant$hpd_97.5 < -rope | quant$hpd_2.5 > rope)
out_rope <- quant$protein[out_rope_log]

report_table <- quant[out_rope_log, c("protein", "mean", "hpd_2.5", "hpd_97.5", "n_peptides")]
colnames(report_table) <- c("Protein", "log2FC", "HPDI start", "HPDI end", "Peptides")

rws <- seq(1, (nrow(report_table)-1), by = 2)
col <- rep("\\rowcolor[gray]{0.95}", length(rws))

xtable(report_table,
       digits=2,
       booktabs = TRUE,
       add.to.row = list(pos = as.list(rws)),
       command=col)


histogram <- ggplot(data = quant, aes(x = mean)) +
  geom_vline(xintercept = c(-rope, rope), linetype="dashed") +
  geom_histogram(bins = 100, position="dodge", fill = "lightgreen", col = "#aaaaaa") +
  labs(x = "log2FC", y = "Count")

histogram
ggsave(filename = file.path(plots_dir, paste0("histogram", opt$suffix, ".png")), plot = histogram)
ggsave(filename = file.path("../../Report/plots/", paste0("histogram", opt$suffix, ".png")), plot = histogram)

dumbbell <- make_dumbbell_plot(quant, split = ", ", ref_val = NULL,n=1000, facet=F)
dumbbell
ggsave(filename = file.path(plots_dir, paste0("dumbbell", opt$suffix, ".png")), plot = dumbbell)
ggsave(filename = file.path("../../Report/plots/", paste0("dumbbell", opt$suffix, ".png")), plot = dumbbell)

width_dist <- ggplot(quant, aes(x = `hpd_97.5` - `hpd_2.5`)) + geom_histogram(fill="Lightgreen", col = "#aaaaaa") + labs(x = "95% HPDI width", y = "Count")
width_dist
ggsave(filename = file.path(plots_dir, paste0("width_dist", opt$suffix, ".png")), plot = width_dist)
ggsave(filename = file.path("../../Report/plots/", paste0("width_dist", opt$suffix, ".png")), plot = width_dist)

width_mean_cor <- ggplot(quant, aes(x = `hpd_97.5` - `hpd_2.5`, y = abs(mean))) + geom_point() +
  labs(y="log2FC", x = "95% HPDI width")
width_mean_cor
cor(quant$hpd_97.5 - quant$hpd_2.5, abs(quant$mean))
ggsave(filename = file.path(plots_dir, paste0("width_mean_cor", opt$suffix, ".png")), plot = width_mean_cor)
ggsave(filename = file.path("../../Report/plots/", paste0("width_mean_cor", opt$suffix, ".png")), plot = width_mean_cor)

bq_comb <- plot_grid(histogram, dumbbell, width_dist, width_mean_cor, nrow=2, labels = "AUTO")
ggsave(filename = file.path(plots_dir, paste0("bq_comb", opt$suffix, ".png")), plot = bq_comb)
ggsave(filename = file.path("../../Report/plots/", paste0("bq_comb", opt$suffix, ".png")), plot = bq_comb, height=8, width = 8)
