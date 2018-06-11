library("dplyr")
library("tidyr")
library("readxl")
library("xtable")
n <- 200
home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
root_dir <- file.path(home_dir, "thesis", "genedata")
input_dir <- file.path(root_dir, "maxlfq/paper")
output_dir <- file.path(root_dir, "maxlfq/quantification")
############################################################################
## Read data from MaxQuant paper
############################################################################

# Read the peptides.txt file
peptides_full <- read.table(file = file.path(input_dir, "peptides.txt"), header = T, sep = "\t", stringsAsFactors = F)
# peptides_full <- peptides_full %>% arrange(Leading.razor.protein)

# Read the supplementary file detailing the Taxonomy - Protein.IDs pairing 
supplementary_full <- read_xlsx(path = file.path(input_dir, "supplementary.xlsx"))
colnames(supplementary_full)[colnames(supplementary_full) == "Protein IDs"] <- "Protein.IDs"
taxonomy <- supplmentary_full[c("Protein.IDs", "Taxonomy")]

# Read the proteinGroups.txt file and keep just the protein groups available in the supplementary file
proteinGroups_full <- read.table(file = file.path(input_dir, "proteinGroups.txt"), header = T, sep = "\t")
proteinGroups <- filter(proteinGroups_full, Protein.IDs %in% supplementary_full$Protein.IDs)
proteinGroups <- inner_join(proteinGroups, supplementary_full, by = "Protein.IDs")
rm(supplementary_full)

# Add the Protein.IDs (protein groups) the peptides belong to.
# This can be done by mapping the Protein.group.IDs column from the peptides.txt file
# to the id column in the proteinGroups.txt!!!! ;)
peptides_full <- merge(peptides_full, select(proteinGroups, id, Protein.IDs),
                       by.x = "Protein.group.IDs",
                       by.y = "id")

# This code will be removed in next release
# peptides_full$Protein.IDs <- NA
# for(i in 1:nrow(peptides_full)) {
#   print(i)
#   mypep <- peptides_full[i,]
#   rowID <- grep(x = proteinGroups$Protein.IDs, pattern = mypep$Leading.razor.protein)
#   rev_rowID <- grep(x = proteinGroups$Protein.IDs, pattern = paste0("REV__", mypep$Leading.razor.protein))
#   rowID <- rowID[!(rowID %in% rev_rowID)]
#   Protein.IDs <- proteinGroups[rowID,"Protein.IDs"]
#   Protein.IDs <- ifelse(length(Protein.IDs) == 1, Protein.IDs, NA)
#   peptides_full[i, "Protein.IDs"] <- Protein.IDs
# }

# peptides_full$Proteins
# peptides <- peptides_full %>% filter(!is.na(Protein.IDs))
peptides <- peptides_full
rm(peptides_full)

# Keep peptides with mapped proteingroup and proteingroups with peptides available
table(peptides$Protein.IDs %in% proteinGroups$Protein.IDs)
table(proteinGroups$Protein.IDs %in% peptides$Protein.IDs)
# proteinGroups[!(proteinGroups$Protein.IDs %in% peptides$Protein.IDs),] %>% View
proteinGroups <- proteinGroups[(proteinGroups$Protein.IDs %in% peptides$Protein.IDs),]
Protein.IDs <- proteinGroups$Protein.IDs
taxonomy <- taxonomy[taxonomy$Protein.IDs %in% Protein.IDs, ]
table(proteinGroups$Protein.IDs %in% peptides$Protein.IDs)

############################################################################
## Extract peptide intensities and sample combinations for this experiment
############################################################################
intensities <- peptides[,grep(x = colnames(peptides), pattern = "Intensity.", value=T)]
colnames(intensities) <- colnames(intensities) %>% gsub(pattern = "Intensity.", replacement = "")


## Impute peptide intensities if it's missing in only one of the replicates but al the other have it
## Maybe we should not be doing this
# for (g in list(1:3, 4:6)) {
#   print(g)
#   imputed_rows <- rowSums(intensities[,g] == 0) == 1
#   observed_values <- intensities[imputed_rows, g] %>% apply(., 1, function(x) sort(x) %>% .[2:3]) %>% t
#   observed_means <- observed_values %>% rowMeans
#   observed_diff <- observed_values[,2] / observed_values[,1]
#   
#   imputed_rows <- which(imputed_rows)[observed_diff < 2]
#   for (row in imputed_rows) {
#     intensities[row, g][,intensities[row, g] == 0] <- mean(intensities[row, g][,intensities[row, g] != 0] %>% unlist) 
#     
#   }
# }

peptides_Protein.IDs <- peptides$Protein.IDs
peptides_sequence <- peptides$Sequence

rowSums(intensities[,1:3] == 0) %>% table
rowSums(intensities[,4:6] == 0) %>% table

combinations <- combn(x = colnames(intensities), m = 2)
combinations_names <- apply(X = combinations, MARGIN = 2, FUN = function(x) paste(x, collapse="/"))

############################################################################################
## Compute the log2 FC based on the spectral counts, the intensities and the LFQ method
############################################################################################

## These are the three methods compared in the MaxLFQ paper
log2ratioLFQ <- log2(rowMeans(proteinGroups[,paste0("LFQ intensity H", 1:3)]) / rowMeans(proteinGroups[,paste0("LFQ intensity L", 1:3)]))
log2ratioIntensities <- log2(rowMeans(proteinGroups[,paste0("Intensity H", 1:3)]) / rowMeans(proteinGroups[,paste0("Intensity L", 1:3)]))
log2ratioSC <- log2(rowMeans(proteinGroups[,paste0("MS/MS Count H", 1:3)]) / rowMeans(proteinGroups[,paste0("MS/MS Count L", 1:3)]))


############################################################################################
## Compute the log10 of the peptide-median summed intensity protein-wise 
############################################################################################

i <- 1
for (protein in Protein.IDs) {
  protein_pos <- which(peptides$Protein.IDs == protein)
  summed_peptide_intensities[i] <- median(intensities[protein_pos, reference_samples] %>% rowSums())
  i <- i + 1
}
log10sumI <- log10(summed_peptide_intensities)

############################################################################################
## Compile dataset
############################################################################################

peptides_data <- peptides %>% select(Amino.acid.before, Sequence, Amino.acid.after,
                    Missed.cleavages, Mass, Protein.IDs, Gene.names, Charges,
                    Intensity:Intensity.L3, PEP, Score)

peptides_data <- inner_join(peptides_data, select(proteinGroups, Protein.IDs, Taxonomy), by = "Protein.IDs")


peptides_data <- peptides_data %>% arrange(Taxonomy, Protein.IDs, Sequence)
peptides_data %>% group_by(Protein.IDs) %>% summarise(count = n()) %>% .$count %>%
  table
peptides_data %>% group_by(Protein.IDs) %>% summarise(count = n()) %>% .$count %>%
  hist(., breaks = 100)


dataset_summary_table <- left_join(peptides_data %>% group_by(Taxonomy) %>% summarise(Peptides = n()),
          peptides_data %>% group_by(Taxonomy) %>% summarise(Proteins = length(unique(Protein.IDs))),
          by = "Taxonomy"
          ) %>% xtable

peptides_data$Protein.IDs %>% unique %>% strsplit(., split = ";") %>% unlist %>% duplicated %>% table
peptides_data$Protein.IDs %>% unique %>% length

############################################################################################
## Export to tsv files
############################################################################################

write.table(x = intensities, file = file.path(output_dir, "maxlfq_peptide_intensities.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = peptides_Protein.IDs, file = file.path(output_dir, "peptides_Protein.IDs.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = peptides_sequence, file = file.path(output_dir, "peptides_sequence.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = combinations, file = file.path(output_dir, "combinations.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = combinations_names, file = file.path(output_dir, "combinations_names.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)

write.table(x = log2ratioLFQ, file = file.path(output_dir, "log2ratioLFQ.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = log2ratioIntensities, file = file.path(output_dir, "log2ratioIntensities.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = log2ratioSC, file = file.path(output_dir, "log2ratioSC.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = log10sumI, file = file.path(output_dir, "log10sumI.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)


write.table(x = Protein.IDs, file = file.path(output_dir, "Protein.IDs.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = taxonomy, file = file.path(output_dir, "Taxonomy.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

write.table(x = peptides_data, file = file.path(output_dir, "peptides_data.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
