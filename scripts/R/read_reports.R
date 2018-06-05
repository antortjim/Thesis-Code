library("dplyr")
library("optparse")
library("ggplot2")

option_list = list(
  make_option(c("--root_dir"), type="character", default="//hest/aoj//thesis/genedata"),
  make_option(c("--exp_name"), type="character", default="thp1"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

input_dir <- ifelse("input_dir" %in% names(opt), opt$input_dir, file.path(root_dir, exp_name, "peptideShaker_out/PSM_reports"))

PSM_report <- read.table(file = file.path(input_dir, "PD7505-GDTHP1-A_C2_Default_PSM_Report.txt"), header = T, sep = "\t", stringsAsFactors = F)
match_report <- read.table(file = file.path(input_dir, "mbr_output", "PD7505-GDTHP1-A_C2_Default_PSM_Report_match.txt"), header = T, sep = "\t", stringsAsFactors = F)

match_report[c(tail(which(match_report$matched == 0)), head(which(match_report$matched == 1))), ] %>% View


protein_report <- read.table(file.path(home, reports_dir), "test_experiment_sample_1_1_Default_Protein_Report_clean.txt"),
                             sep = "\t", header=T, na.strings = "")
my_report <- read.table(file.path(home, reports_dir), "my_report_clean.txt"),
                        sep = "\t", header=T, na.strings = "")
peptide_report <- read.table(file.path(home, reports_dir), "test_experiment_sample_1_1_Default_Peptide_Report_clean.txt"),
                             sep = "\t", header=T, na.strings = "")

# peptide_report <- read.table(file.path(home, reports_dir), "peptide_report_clean.txt"),
#                              sep = "\t", header=T, na.strings = "")

psm_report <- read.table(file.path(home, reports_dir), "test_experiment_sample_1_1_Default_PSM_Report_clean.txt"),
                         sep = "\t", header=T, na.strings = "", stringsAsFactors = F)

mz_i <- read.table(file = file.path(home, reports_dir), "thesis/genedata/thp1/data/mgf/mz_intensity.txt"))
scan <- read.table(file = file.path(home, reports_dir), "thesis/genedata/thp1/data/mgf/scans.txt"))
precursor <- data.frame(Spectrum.Scan.Number = scan[,1], mz = mz_i[,1], intensity = mz_i[,2])

psm_report <- left_join(psm_report, select(precursor, Spectrum.Scan.Number, intensity), by = "Spectrum.Scan.Number")

new_psm_report <- psm_report[NULL,]
for (i in 1:nrow(psm_report)) {
  print(i)
  proteins <- strsplit(psm_report[i,"Protein.s."] %>% as.character(), split = ", ") %>% unlist
  current <- psm_report[rep(i, times = length(proteins)),]
  current$Protein.s. <- proteins
  new_psm_report <- rbind(new_psm_report, current)
}

unique_psms <- (new_psm_report %>% group_by(Spectrum.Title) %>% summarise(count = n()) %>% filter(count == 1) %>% .$Spectrum.Title)
unique_psm_report <- new_psm_report %>% filter(Spectrum.Title %in% unique_psms) 
new_psm_report %>% filter(Spectrum.Title %in% unique_psms) %>% .$Protein.s. %>% table %>% sort

si <- unique_psm_report %>% group_by(Protein.s., Sequence) %>% summarise(intensity = sum(intensity)) %>%
  group_by(Protein.s.) %>% summarise(intensity = sum(intensity))

si2 <- si$intensity
names(si2) <- si$Protein.s.
si <- si2

mean_si <- mean(si)

si_actin <- si["O43707"]

si_act <- si / si_actin

si_mpi <- si / mean_si

si_gi <- si / sum(si)

si_gi_df <- data.frame(name = names(si_gi), si_gi = si_gi)

compare <- left_join(select(data, name, amount), si_gi_df, by = "name")
plot(compare$amount, compare$si_gi)

linear_model <- lm(formula = amount ~si_gi, data = compare)

summary_linear_model <- summary(linear_model)
summary_linear_model$adj.r.squared

library("seqinr")
database <- read.fasta(file = "thesis/genedata/thp1/databases/all_concatenated_target_decoy.fasta", seqtype = "AA")
annotations <- getAnnot(database)
sequences <- getSequence(database)
protein_ids <- annotations %>% lapply(function(x) strsplit(x, split = "\\|") %>% unlist %>% .[2]) %>% unlist
protein_lengths <- sequences %>% lapply(length) %>% unlist
names(protein_lengths) <- protein_ids

si_n <- si_gi / protein_lengths[names(si_gi)]

setwd("../..")



peptide_shaker <- my_report %>% select(name = Main.Accession, nsaf = Spectrum.Counting.NSAF..fmol.)

compare_data <- full_join(peptide_shaker, data[data$file == 7, c("name", "amount")], by = "name")


compare_data[,-1] %>% apply(X = ., MARGIN = 1, FUN = function(x) sum(which(is.na(x)))) %>% unlist %>% table
compare_data <- compare_data[complete.cases(compare_data),]


plot(compare_data$nsaf, compare_data$amount, pch=16)
abline(a = 0,b =1)
