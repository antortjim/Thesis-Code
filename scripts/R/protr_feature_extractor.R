setwd("thesis/genedata/maxlfq/paper")
peptides_data <- read.table(file = "peptides_data.tsv", header=T, sep = "\t", stringsAsFactors = F)
library("protr")

peptides_sequence <- peptides_data$Sequence
length(peptides_sequence)

# Ensure standard 20 AA code
peptides_sequence <- peptides_sequence[sapply(peptides_sequence, protcheck)]
length(peptides_sequence)

peptides_length <- sapply(peptides_sequence, nchar)

hist(peptides_length)
min_length <- min(peptides_length) - 1

# This step takes a lot of time!!! ~ 30 mins
system.time(ApseAAC <- t(sapply(peptides_sequence, function(x) extractAPAAC(x = x, lambda = min_length))))
