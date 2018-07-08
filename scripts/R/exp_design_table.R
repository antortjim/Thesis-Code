library("dplyr")
library("optparse")
library("ggplot2")
library("readxl")
library("stringr")
library(xtable)

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq")
  # make_option(c("--input_dir"), type="character"),
  # make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
input_dir <- opt$input_dir
output_dir <- opt$output_dir

experimental_design <- read.table(file = file.path(root_dir, exp_name, "data", "experimental_design.tsv"), sep = "\t", header=T, stringsAsFactors = F)
experimental_design <- experimental_design %>% mutate(Replicate = factor(Replicate))
experimental_design <- rbind(experimental_design, data.frame(Name=NA, Group="6_H", Fraction=6, Experiment="H", Replicate=1))

exp_table <- select(experimental_design, Experiment, Replicate, Fraction) %>%
  arrange(Experiment, Fraction, Replicate) %>% filter(Replicate == 1) %>% select(-Replicate)

exp_table <- cbind(
  exp_table %>% filter(Experiment == "H"),
  exp_table %>% filter(Experiment == "L")
)


xtable::xtable(exp_table)
