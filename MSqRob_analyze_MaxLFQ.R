###############################
## Load libraries
###############################
library(MSqRob)
library(limma)
library(Biobase)
library(dplyr)
library("here")
library("Biostrings")
library(protr)
library(tidyr)
library(norm)
library(readxl)
library("ggplot2")

aacode_number <- 1:length(Biostrings::AA_ALPHABET)
names(aacode_number) <- Biostrings::AA_ALPHABET

#################################
## Load data
#################################
setwd(file.path(here(), "MEGA/Master/Thesis/Code/"))
file_peptides_txt <- "MaxLFQ/peptides.txt"
file_proteinGroups <- "MaxLFQ/proteinGroups.txt"
## Read the peptides.txt file
peptidesMaxLFQ <- import2MSnSet(file_peptides_txt, filetype="MaxQuant", remove_pattern=TRUE)
exp_annotation <- read.table(file = "exp_annotation.tsv", header = T, sep = "\t")
###################################


# Standard preprocessing of MaxQuant data
peptidesMaxLFQ_preprocessed <- preprocess_MaxQuant(peptidesMaxLFQ, # the MSnset containing the peptides results
                                      accession="Proteins", # column containing the protein identifiers, required if smallestUniqueGroups is TRUE
                                      exp_annotation=exp_annotation, # experimental annotation to take batch effect into account
                                      logtransform=TRUE, # log transform intensities (to make them more symmetrical  https://www.sciencedirect.com/science/article/pii/S1874391917301239?via%3Dihub)
                                      base=2, # base to use in the log transform
                                      # normalisation="quantiles", # remove systematic bias across samples (runs)
                                      normalisation="none",
                                      smallestUniqueGroups=TRUE, # occam razor for the protein groups
                                      useful_properties=c("Proteins","Sequence"), # columns to keep in further analysis
                                      filter=c("Potential.contaminant","Reverse"), # remove entries positive for these fields
                                      remove_only_site=TRUE, # true if we want to discard proteins identified with modified peptides
                                      file_proteinGroups=file_proteinGroups, # this needs to be passed if remove_only_site is True
                                      filter_symbol="+", # if +, the field is positive
                                      minIdentified=2 # the peptide needs to have been identified at least twice in the dataset
                                      )

# Introduce into the MSqRob pipeline
# MSnSet -> protdata
system.time(proteinsMaxLFQ <- MSnSet2protdata(peptidesMaxLFQ_preprocessed, accession="Proteins"))

data <- getData(proteinsMaxLFQ)
annotations <- getAccessions(proteinsMaxLFQ)
# data_df <- lapply(1:length(data), function(i) {print(i); cbind(data[[i]], protein = annotations[i])}) %>% do.call(rbind, .)
# write.table(x = data_df, file = "data/peptides_df.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


############################################################################################################################
### Avanced and optional: check which models would fail to converge by building a linear mixed effects model.
############################################################################################################################

# MSqRob always tries to fit a model, but some models are overparameterized
# because too many parameters are fit on too few observations.
# These models have convergence problems and can be removed from the data prior to estimating p-values.
# This is only relevant for models that perform shrinkage or use any kind of random effect.

convergence <- lapply(1:length(data), function(i){print(i); return(tryCatch(lme4::lFormula(
  formula("quant_value~1+(1|condition)+(1|Sequence)"), data[[i]], control = lme4::lmerControl()
), error=function(e){
  return(NA)
}))})

na_indices <- which(is.na(convergence))
na_indices %>% length
data_converge <- data[-na_indices]

# # Filtered peptides: a character with the sequences that we wanna use downstream
#   # Take those available in the converged object
# filtered_peptides <- lapply(data_converge, function(x) x$Sequence) %>% unlist %>% unique %>% as.character
# annotations_converge <- as.character(annotations)[-na_indices]
# proteinsMaxLFQ_converge <- proteinsMaxLFQ[-na_indices]
# data_annotated <- lapply(1:length(data_converge), function(i) {
#   print(i);
#   cbind(data_converge[[i]], protein = annotations_converge[i])
#   }) %>% do.call(rbind, .)
# proteins <- data_annotated$protein %>% as.character %>% unique
# Export a list of the available proteins
# write.table(file = "proteins.txt", x = proteins, quote = F, row.names = F, col.names = F)

###############################################################################
### Add annotations for each protein group
### The annotations_df.tsv file will have a column called protein with protein ids
### similar to those in data_annotated, but with extra columns bringing more information
### about this protein ids
### No missing data is present in annotations_df, thus if not all the annotations
### could be found for a protein group, it will be dropped
### This behaviour can be changed by changing the add_annotations.R file
### which currently runs complete.cases() on the output before exporting to a tsv.
###############################################################################
# source("check_organism.R")
# proteins <- unique(data_annotated$protein) %>% as.character
# organism_check <- check_organism(proteins)
# Taxonomy <- data.frame(protein = proteins, taxonomy = organism_check[[2]])

# data_annotated_clean <- inner_join(data_annotated, Taxonomy, by = "protein")
# # data_annotated %>% group_by(protein) %>% summarise(missing = all(is.na(taxon))) %>% .$missing %>% table
# dropped_peptides <- data_annotated[!(as.character(data_annotated$Sequence) %in% as.character(data_annotated_clean$Sequence)), "Sequence"]
# filtered_peptides <- filtered_peptides[!(filtered_peptides %in% dropped_peptides)]

# ## Extract sequence features
# #############################################################################################################################
# ## Extract the sequence from the featureData@data slot
# sequences <- peptidesMaxLFQ@featureData@data %>%
#   # Select the actual sequence pluts the aminoacid before, after and the number of missed cleavages
#   select(Amino.acid.before, Sequence, Amino.acid.after, Missed.cleavages) %>%
#   # Keep those peptides that were filtered upstream (i.e belonging to proteins that converge and proteins for which full annotation is available)
#   filter(Sequence %in% filtered_peptides)
# sequences$Sequence <- sequences$Sequence %>% as.character
# write.table(x = sequences, file = "sequences.txt", row.names = F, col.names=T, quote=F, sep = "\t")
# write.table(x = data_annotated_clean, file = "data_annotated_clean.tsv", row.names = F, col.names=T, quote=F, sep = "\t")
# #### CALL extract_features.R
# #############################################################################################################################

# takes a few: 401 seconds

# empirical bayes estimation of variance makes the protein-wise-estimated error term variance more robust.
# this is needed when few peptides are sampled from a protein, which could lead to inflated underlying variance estimates
# i.e the sample variance is much bigger than the actual variance
# it is enabled by passing squeezeVar=TRUE

# fixed variables are those where all possible levels are included in the experiment
# random are those where not all are included, probably because there's infinite amount.
# the levels included in random variables are thus "considered to be drawn at random from a broader near-infinite population"
#

# weights = Huber/ NULL: if NULL, the normal ridge regression loss is used i.e the squared error loss + the ridge regression terms
# if Huber, the loss used differs for outlier datapoints, defined as points standing from the average trend of the dataset as more as delta
# https://en.wikipedia.org/wiki/Huber_loss

# shrinkage.fixed indicates which fixed effects should be shrunken using ridge regression (whose function is to penalize high values of parameters
# to make them less sensitive to slight changes in data)
# if NULL, the default, all fixed effects except the intercept are shrunken
# otherwise, it's an integer vector of 1/0 and length equal to fixed, where the ith element states
# if the ith fixed effect should be shrunken (1) or not (0).

# takes a long time -> 1000 seconds (quarter of an hour)
system.time(modelMaxLFQ_RR <- fit.model(protdata=proteinsMaxLFQ_converge,
                                       response="quant_value",
                                       fixed=c("condition"),
                                       random=c("Sequence","run"),
                                       shrinkage.fixed=NULL,
                                       weights="Huber",
                                       squeezeVar=TRUE # empirical bayes estimation of protein variance
                                       ))



# attr(getModels(modelMaxLFQ_RR[2]),"MSqRob_levels")
# paste0("condition", exp_annotation$condition)

L <- makeContrast(contrasts=c("conditionH-conditionL"),
                  levels=c("conditionH","conditionL"))

# with the fitted model, perform the contrast testing
RSqM <- test.protLMcontrast(modelMaxLFQ_RR, L)
RSqM_old <- RSqM
source("check_organism.R")
proteins <- rownames(RSqM)
organism_check <- check_organism(proteins)
taxon <- organism_check[[2]]
RSqM_mapping <- data.frame(protein = proteins, taxonomy = taxon)
RSqM$protein <- rownames(RSqM)
RSqM <- inner_join(RSqM, RSqM_mapping, by = "protein")
rownames(RSqM) <- NULL


# adjust the pvalues using Benjamini Hochberg
RSqM_adjust <- prot.p.adjust(RSqM)
RSqM_signif <- prot.signif(RSqM_adjust)

RSqM_signif_clean <- RSqM_signif %>% filter(!is.na(taxonomy))

ggplot(RSqM_signif_clean, aes(x = estimate, fill = taxonomy)) + geom_histogram(position = "dodge", bins = 60)
ggsave("plots/MSqRob_quantification_maxquant_data.png")
ggplot(RSqM_signif_clean, aes(x = estimate)) + geom_histogram(position = "dodge")
write.table(x = RSqM_signif, file = paste0(file.path("data", "RSqM_signif", "_", file_peptides_txt)), quote = F,row.names = F, col.names = T, sep = "\t")