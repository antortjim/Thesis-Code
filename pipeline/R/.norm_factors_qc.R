library("dplyr")
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

norm_factors <- read.table(file.path(root_dir, exp_name, "peptideShaker_out", "PSM_Reports", "output_moff_RAW", "normalization_factors.tsv"))[,1]
sample_names <- read.table(file.path(root_dir, exp_name, "peptideShaker_out", "PSM_Reports", "output_moff_RAW", "peptide_summary_intensity_moFF_run.tab"),
                           nrows = 1, header = F)[1,][-(1:2)] %>% unlist %>% as.character %>%
  gsub(pattern = "sumIntensity_", replacement = "")

experimental_design <- read.table(file.path(root_dir, exp_name, "data", "experimental_design.tsv"), header=T)
norm_factors <- data.frame(norm_factor = norm_factors, Name = sample_names)

norm_factors <- full_join(experimental_design, norm_factors,by="Name")

library("ggplot2")

ggplot(data = norm_factors, aes(x = norm_factor, fill = Experiment)) + geom_histogram(position="dodge")

ggplot(data = norm_factors, aes(x = norm_factor)) +
  geom_histogram(position="dodge", bins=60)
