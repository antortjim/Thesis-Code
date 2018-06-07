.libPaths(new="/nfs/home/aoj/.R")
print(.libPaths())
library("optparse")
library("dplyr")

option_list = list(
  make_option(c("--root_dir"), type="character", default="//hest/aoj//thesis/genedata"),
  make_option(c("--exp_name"), type="character", default="thp1_test"),
  make_option(c("--sample_names"), type="character", default="PD7505-GDTHP1-A_C2|PD7506-GDTHP1-A_C1|PD7507-GDTHP1-A_C2")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
sample_names <- opt$sample_names %>% strsplit(split = "\\|") %>% unlist


# if(interactive()) {
#   root_dir <- "//hest/aoj//thesis/genedata"
#   exp_name <- "thp1_test"
#   sample_names <- c("PD7505-GDTHP1-A_C2","PD7506-GDTHP1-A_C1","PD7507-GDTHP1-A_C2")
# } else {
#   arguments <- commandArgs(trailingOnly=TRUE)
# 
#   root_dir <- arguments[[1]]
#   exp_name <- arguments[[2]]
#   sample_name <- arguments[[3]] %>% strsplit(split = "\\|") %>% unlist
# 
# }


for (sample_name in sample_names) {
  print(sample_name)
  custom_psm_report_path  <- paste0(root_dir, "/", exp_name, "/peptideShaker_out/custom_reports/", sample_name, "_Custom_PSM_Report.txt")
  default_psm_report_path  <- paste0(root_dir, "/", exp_name, "/peptideShaker_out/PSM_reports/", sample_name, ".txt")
  matched_psm_report_path  <- paste0(root_dir, "/", exp_name, "/peptideShaker_out/PSM_reports/mbr_output/", sample_name, "_match.txt")
  source_psm_report_path <- default_psm_report_path
  
  
  con <- file(source_psm_report_path,"r")
  old_header <- readLines(con,n=1)
  close(con)
  custom_header <- paste0(old_header, "\tms1_intensity")
  
  source_psm_report <- read.table(file = source_psm_report_path, header = T,sep = "\t", comment.char = "")
  ms1 <- read.table(file = paste0(root_dir, "/", exp_name, "/data/ms1/", sample_name, "_ms1_table.tsv"), header = T,sep = "\t", comment.char = "")
  colnames(ms1)[1] <- "spectrum.scan.number"
  colnames(source_psm_report) <- colnames(source_psm_report) %>% tolower
  
  custom_psm_report <- left_join(source_psm_report, select(ms1, spectrum.scan.number, pepmass_intensity), by = "spectrum.scan.number")
  
  write.table(x = custom_header, file = custom_psm_report_path, col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(x = custom_psm_report, file = custom_psm_report_path, col.names = F, append = T, row.names = F, quote = F, sep = "\t")
}
