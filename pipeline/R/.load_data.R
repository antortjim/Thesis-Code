library("dplyr")
library("optparse")
library("ggplot2")
library("readxl")
library("stringr")
library(xtable)

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="thp1")
  # make_option(c("--input_dir"), type="character"),
  # make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
input_dir <- opt$input_dir
output_dir <- opt$output_dir

read_ms <- function(filename) {
  data_dir <- file.path(root_dir, exp_name, "data", "xlsx")
  # Reads the excel file in filename
  file_path <- file.path(data_dir, filename)
  print(filename)
  dat <- read_excel(file_path)
  metadata_fields <- c("Source Type", "Group", "Description", "Maximum RT Shift", "Flags", "Search Results", "Enzyme", "Fraction")
  metadata_rows <- which(dat$Name %in% metadata_fields)
  metadata <- dat[metadata_rows, 1:2]
  dat <- dat[-metadata_rows,]
  
  dat <- dat[-(dat$Name %>% grep(pattern = "con_| : ", x = .)), ]
  
  colnames(dat)[2] <- "Amount"
  colnames(dat) <- colnames(dat) %>% make.names() %>% tolower
  numeric_cols <- apply(dat, 2, function(x) as.numeric(x) %>% is.na %>% all %>% !.) %>% which
  dat[,numeric_cols] <- apply(dat[,numeric_cols], 2, as.numeric)
  # Add file id: file i means files[i]
  dat$Name <- str_match(string = filename, pattern = "\\w \\[(.*)\\]\\.xlsx")[,2]
  dat <- dat[,c("name", "amount", "Name", "protein.p.value")]
  colnames(dat)[colnames(dat) == "name"] <- "Protein.ID"
  colnames(dat)[colnames(dat) == "amount"] <- "Quantification"
  colnames(dat)[colnames(dat) == "protein.p.value"] <- "PEP"
  
  return(dat)
}



data_list <- list()
files <- list.files(path = file.path(root_dir, exp_name, "data", "xlsx"), "*.xlsx")

experimental_design <- read.table(file = file.path(root_dir, exp_name, "data", "experimental_design.tsv"), sep = "\t", header=T, stringsAsFactors = F)
experimental_design <- experimental_design %>% mutate(Replicate = factor(Replicate))

i <- 1
for (filename in files) {
  data_list[[i]] <- read_ms(filename)
  i <- i + 1
}

data <- do.call(rbind, data_list)
data <- data %>% full_join(experimental_design, by = "Name")


write.table(x = data, file = file.path(root_dir, exp_name, "data", "genedata_processed", "dataset.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = data$Protein.ID %>% unique, file = file.path(root_dir, exp_name, "data", "genedata_processed", "ids.txt"), quote = F, row.names = F, col.names=F)
