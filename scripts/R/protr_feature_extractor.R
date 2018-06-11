setwd("thesis/genedata/maxlfq/paper")
peptides_data <- read.table(file = "peptides_data.tsv", header=T, sep = "\t", stringsAsFactors = F)
library("protr")
library("dplyr")

peptides_sequence <- peptides_data$Sequence
length(peptides_sequence)

# Ensure standard 20 AA code
peptides_sequence <- peptides_sequence[sapply(peptides_sequence, protcheck)]
length(peptides_sequence)

peptides_length <- sapply(peptides_sequence, nchar)

hist(peptides_length)
min_length <- min(peptides_length) - 1

# This step takes a lot of time!!! ~ 30 mins
aminoacids <- peptides_sequence %>% paste(., collapse="") %>% strsplit(., split = "") %>% table

system.time(ApseAAC <- t(sapply(peptides_sequence, function(x) extractAPAAC(x = x, lambda = min_length))))
write.table(x = ApseAAC, file = "apseaac_features.tsv", sep = "\t", quote = F, row.names = F, col.names = F)
