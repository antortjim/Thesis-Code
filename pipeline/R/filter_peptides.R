library("dplyr")
library("tidyr")

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj",
                   ifelse(Sys.info()["user"] == "aoj", "/z/home/aoj",
                          "/home/antortjim/"))


data_dir <- ifelse(Sys.info()["user"] == "aoj", "thesis/genedata/maxlfq/",
                   "MEGA/Master/Thesis/Code/scripts/data")


peptide_summary <- read.table(file.path(home_dir, data_dir,
                                        "output_moFF_RAW/peptide_summary_intensity_moFF_run.tab"
                                        ),
                              header=T, sep = "\t")

peptide_summary_log2 <- cbind(peptide_summary[,1:2], log2(peptide_summary[, -(1:2)]))

colnames(peptide_summary_log2) <- colnames(peptide_summary_log2) %>% gsub(pattern = "sumIntensity_", replacement = "")

exp_annot <- read.table(file.path(home_dir, data_dir, "experimental_design.tsv"), header = T, stringsAsFactors = F)

groups <- exp_annot$Group %>% unique %>% as.character

group_check <- matrix(ncol=length(groups), nrow=nrow(peptide_summary))
for (i in 1:length(groups)) {
  print(i)
  current_fraction <- peptide_summary_log2[,exp_annot %>% filter(Group == groups[i]) %>% .$Name]
  cvs <- apply(current_fraction,1,function(x) sd(x)/mean(x))
  cvs[is.na(cvs)] <- 0
  group_check[,i] <- cvs < 1
}

peptide_summary_log2_clean <- peptide_summary_log2[(group_check %>% apply(., 1, all)),]


peptide_summary_log2_clean[,exp_annot %>% filter(Group == groups[i]) %>% .$Name] %>% View
