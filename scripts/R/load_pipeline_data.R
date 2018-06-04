library("dplyr")
library("tidyr")
library("readxl")
library("here")
library("xtable")
library("optparse")
library("stringr")

option_list = list(
  make_option(c("--root_dir"), type="character", default="//hest/aoj//thesis/genedata"),
  make_option(c("--exp_name"), type="character", default="thp1_test")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

custom_reports_dir <- file.path(root_dir, exp_name, "peptideShaker_out/custom_reports")
output_dir <- file.path(root_dir, exp_name, "quantification")

PSM_reports <- list.files(custom_reports_dir)
intensities <- data.frame(Sequence = "")
peptides <- data.frame()

i <- 1
sample_names <- str_match(string = PSM_reports, pattern = "(.*)_Custom_PSM_Report.txt")[,2]

for (PSM_rep in PSM_reports) {
  sample_name <- sample_names[i]
  
  pept <- read.table(file = file.path(custom_reports_dir, PSM_rep),
             sep = "\t", header=T, comment.char = "", stringsAsFactors = F)
  
  peptides_non_contaminant <- pept[-(grep(pattern = "_CONTAMINANT", pept$Protein.s)),] %>%
    select(Protein.s.:Modified.Sequence,
           #Spectrum.Scan.Number, RT,
           Measured.Charge, ms1_intensity)
  
  peptides <- rbind(peptides,
                    cbind(peptides_non_contaminant, sample_name = sample_name)
  )
  
  psm_rep <- peptides_non_contaminant %>%
      select(Sequence, Modified.Sequence, Measured.Charge, ms1_intensity) %>%
    group_by(Sequence) %>% summarise(summed_intensity = sum(ms1_intensity))
  
  colnames(psm_rep)[colnames(psm_rep) == "summed_intensity"] <- sample_name
  # psm_rep_list[[i]] <- psm_rep
  
  i <- i + 1
}


pp <- peptides %>% select(-Modified.Sequence) %>%
  unite(temp, Protein.s.:Measured.Charge)

pp$temp %>% duplicated
pp <- pp %>% spread(sample_name, ms1_intensity)

temp <- pp$temp %>% strsplit(x = ., split = "_") %>% do.call(what = rbind, args = .)

colnames(temp) <- peptides %>% select(-Modified.Sequence) %>% select(Protein.s.:Measured.Charge) %>% colnames
pp <- cbind(temp, select(pp, -temp))
rm(temp)
peptides <- pp
colnames(peptides)[colnames(peptides) == "Protein.s."] <- "Protein.IDs"
Protein.IDs <- peptides$Protein.IDs %>% as.character

intensities <- peptides[,sample_names]
combinations <- combn(x = colnames(intensities), m = 2)
combinations_names <- apply(X = combinations, MARGIN = 2, FUN = function(x) paste(x, collapse="/"))


############################################################################################
## Export to tsv files
############################################################################################

write.table(x = intensities, file = file.path(output_dir, "intensities.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = combinations, file = file.path(output_dir, "combinations.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = combinations_names, file = file.path(output_dir, "combinations_names.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = peptides, file = file.path(output_dir, "peptides.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = Protein.IDs, file = file.path(output_dir, "Protein.IDs.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
