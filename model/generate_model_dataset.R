library(GenomicRanges)
library(seqinr)
library(dplyr)
library(readxl)
library(Biostrings)
library(MSqRob)
library(MSnbase)
library(tidyr)
library(optparse)

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--data_dir"), type="character", default="data"),
  make_option(c("--extract-features"), action="store_true", default=FALSE,
              help="Extract sequence features")
  # right now it breaks if extract-features = FALSE
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data_dir <- opt$data_dir
extract_features <- opt$extract_features

source("check_organism.R")
strsplit <- base::strsplit

file_peptides_txt <- "data/peptides.txt"

# Import MSnSet
peptidesMaxLFQ <- import2MSnSet(file_peptides_txt, filetype="MaxQuant", remove_pattern=TRUE)

file_proteinGroups <- "data/proteinGroups.txt"
exp_annotation <- read.table(file = "data/exp_annotation.tsv", header = T, sep = "\t")
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
data_df <- data %>% do.call(rbind, .)
accessions <- getAccessions(proteinsMaxLFQ)
data_df$Proteins <- rep(accessions, lapply(data, nrow) %>% unlist)
peptides <- data_df


## SEQUENCE FEATURE EXTRACTION
#############################################################################################

if(extract_features) {
  ecoli_fasta <- read.fasta(file = "fasta/ecoli.fasta",seqtype = "AA")
  human_fasta <- read.fasta(file = "fasta/homo_sapiens.fasta",seqtype = "AA")
  # peptides_head <- head(peptides, 10)
  
  peptides <- peptides %>% select(Sequence, Proteins)
  peptides <- peptides[(!peptides %>% duplicated),]
  
  organism_check <- check_organism(proteins = peptides$Proteins %>% as.character() %>% unique, split = ";")
  peptides <- left_join(data.frame(Proteins = unique(peptides$Proteins), taxon = organism_check[[2]]) %>% filter(!is.na(taxon)), peptides, by = "Proteins")
  
  ecoli_prots <- getAnnot(ecoli_fasta) %>%
    lapply(., function(x) strsplit(x, split = "\\|") %>% unlist %>% .[2])
  human_prots <- getAnnot(human_fasta) %>%
    lapply(., function(x) strsplit(x, split = "\\|") %>% unlist %>% .[2])
  
  
  peptides %>% group_by(taxon) %>% summarise(count=n())
  peptides %>% group_by(taxon) %>% summarise(count=length(unique(Proteins)))
  peptides <- peptides %>% arrange(taxon, Proteins)
  peptides$Sequence <- as.character(peptides$Sequence)
  peptides$Proteins <- as.character(peptides$Proteins)
  
  
  proteins_list <- list("Homo sapiens" = human_prots, "Escherichia coli (strain K12)" = ecoli_prots)
  fastas_list <- list("Homo sapiens" = human_fasta, "Escherichia coli (strain K12)" = ecoli_fasta)
  
  
  match_list <- list()
  not_found_list <- list()
  
  # i <- 1
  n_prots <- 1
  k <- 1
  organism <- "Homo sapiens"
  peptides$Proteins %>% unique %>% strsplit(., split=";") %>% unlist %>% table %>% sort %>% table
  peptides$Proteins %>% unique %>% length
  peptides$prot <- peptides$Proteins %>% strsplit(., split=";") %>% lapply(., function(x) x[1]) %>%
    lapply(., function(x) strsplit(x, split = "-") %>% unlist %>% .[1]) %>%
    unlist
  
  # Takes a couple of minutes
  ################################################################################
  for (organism in unique(peptides$taxon)) {
    taxon_peptides <- peptides %>% filter(taxon == organism)
    taxon_proteins <- proteins_list[[organism]]
    taxon_fasta <- fastas_list[[organism]]
    proteins <- taxon_peptides$prot %>% unique
    for (p in proteins) {
      p_peptides <- taxon_peptides %>% filter(prot == p)
      # split the protein group in the constituent protein ids and remove isoforms
      # ps <- strsplit(protein_group, split=";") %>% unlist %>% lapply(., function(x) strsplit(x, split = "-") %>% unlist %>% .[1]) %>% unlist %>% unique
      fasta_pos <- which(taxon_proteins %in% p)
      if(length(fasta_pos) > 0) {
        protein_sequence <- getSequence(taxon_fasta[fasta_pos])[[1]] %>%
          paste(., collapse="") %>%
          # lapply(., function(x) paste(x, collapse="")) %>% unlist
          AAString
        # Biostrings::(x = protein_sequences)
        # for (pep in p_peptides$Sequence) {
        # print(pep)
        match_range_list <- lapply(p_peptides$Sequence, function(pep) matchPattern(pattern = pep, subject = protein_sequence))
        match_list[[p]] <- match_range_list
        
        # iranges_list[[i]] <- match_range
        # start(match_range) <- max(1, start(match_range) - 15)
        # end(match_range) <- min(width(protein_sequences[1]), end(match_range) + 15)
        # window_15[[pep]] <- protein_sequences[[1]][match_range]
        # i <- i + 1
        # } # end of peptide
      } else {
        print(paste0(p, " not found"))
        not_found_list[[k]] <- p
        match_list[[p]] <- NULL
        
        k <- k + 1
      }
      if (n_prots %% 500 == 0) print(paste0(n_prots, " analyzed"))
      n_prots <- n_prots + 1
    } # end of protein
  } # end of taxon
  ################################################################################
  match_list %>% length
  not_found_list %>% length
  
  # names(match_list)[1]
  (match_list)[[1]][[1]]@subject %>% length
  
  n_rows <- match_list %>% lapply(length) %>% unlist %>% sum
  result <- matrix(nrow=n_rows, ncol=3)
  colnames(result) <- c("pep", "neigh", "prot")
  window_length <- 15
  i <- 1
  k <- 1
  match_list %>% length
  
  # Takes around 5 mins
  ################################################################################
  for(p in names(match_list)) {
    for(m in match_list[[p]]) {
      if (length(m) == 1) {
        pep <- as.character(m)
        start(m) <- max(1, start(m) - window_length);
        end(m) <- min(m@subject %>% length , end(m) + window_length)
        result[i, 1] <- pep
        result[i, 2] <- as.character(m)
        result[i, 3] <- p
        i <- i + 1
      } else {
        print(paste0("Repetitions found for ", as.character(m)[1]))
      }
    }
    k <- k + 1
    if(k %% 100==0) print(k)
  }
  
  result_df <- result %>% as.data.frame()
  head(result_df)

## END
#############################################################################################
# Keep the peptides for which features could be extracted
data_df_subset <- data_df %>% filter(Sequence %in% result[, "pep"])
} else {
  data_df_subset <- data_df
}

# Get organism link
organism_check <- check_organism(proteins = as.character(unique(data_df_subset$Proteins)), split = ";")
taxonomy <- data.frame(Proteins = as.character(unique(data_df_subset$Proteins)), taxon = organism_check[[2]])
data_df_subset <- left_join(data_df_subset, taxonomy, by = "Proteins")

# Make data wide
data_df_subset_wide <- data_df_subset %>% select(-condition) %>% spread(run, quant_value)

# Remove missing data
data_df_subset_wide <- data_df_subset_wide[complete.cases(data_df_subset_wide),]

result_df <- result_df[-(result_df %>% is.na %>% which(arr.ind = T) %>% .[, 1]),]

# nrow(result_df)
nrow(data_df_subset_wide)
colnames(data_df_subset_wide)[1:3] <- c("peptide", "protein", "taxon")
data_df_subset_wide %>% group_by(taxon) %>% summarise(count = n())

## Keep proteins with at least 2 peptides available
solid_proteins <- which((data_df_subset_wide$protein %>% table) >= 2) %>% names
data_df_subset_wide_solid <- data_df_subset_wide %>% filter(protein %in% solid_proteins)


result_df <- result_df[result_df[,"pep"] %in% data_df_subset_wide_solid$peptide,]

# Sort both tables in the same order!!
data_df_subset_wide_solid$peptide <- as.character(data_df_subset_wide_solid$peptide)
data_df_subset_wide_solid <- data_df_subset_wide_solid %>% arrange(peptide)
result_df <- result_df %>% arrange(pep)

all((data_df_subset_wide_solid$peptide == as.character(result_df$pep)))
print(paste0("Exporting data for ", nrow(data_df_subset_wide_solid), " peptides"))


write.table(file = "data/data.tsv", sep = "\t", x = data_df_subset_wide_solid %>% select(-peptide), col.names=T, row.names=F, quote=F)
write.table(file = "data/peptides_neighborhood.tsv", sep = "\t", x = result_df, col.names=T, row.names=F, quote=F)
# save.image("~/MEGA/Master/Thesis/MaxLFQ/extract_advanced_features.RData")
