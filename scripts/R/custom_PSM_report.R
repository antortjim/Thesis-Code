library("dplyr")

arguments <- commandArgs(trailingOnly=TRUE)

root_dir <- arguments[[1]]
exp_name <- arguments[[2]]
sample_name <- arguments[[3]]


default_psm_report_path <- paste0(root_dir, "/", exp_name, "/peptideShaker_out/reports/", sample_name, "_Default_PSM_Report.txt")
custom_psm_report_path  <- paste0(root_dir, "/", exp_name, "/peptideShaker_out/custom_reports/", sample_name, "_Custom_PSM_Report.txt")
new_psm_report_path  <- paste0(root_dir, "/", exp_name, "/peptideShaker_out/PSM_reports/", sample_name, "_Default_PSM_Report.txt")
con <- file(default_psm_report_path,"r")
custom_header <- readLines(con,n=1)
close(con)
custom_header <- paste0(custom_header, "\tms1_intensity")

default_psm_report <- read.table(file = default_psm_report_path, header = T,sep = "\t", comment.char = "")
ms1 <- read.table(file = paste0(root_dir, "/", exp_name, "/data/ms1/", sample_name, "_ms1_table.tsv"), header = T,sep = "\t", comment.char = "")

custom_psm_report <- left_join(default_psm_report, select(ms1, Spectrum.Scan.Number, pepmass_intensity), by = "Spectrum.Scan.Number")

write.table(x = custom_header, file = custom_psm_report_path, col.names = F, row.names = F, quote = F, sep = "\t")
write.table(x = custom_psm_report, file = custom_psm_report_path, col.names = F, append = T, row.names = F, quote = F, sep = "\t")
