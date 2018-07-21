#!/usr/bin/Rscript

library(dplyr)
library(optparse)

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--pepf"), type="character"),
  make_option(c("--protein_file"), type="character", default=""),
  make_option(c("--annotation_df"), type="character"),
  make_option(c("--exp"), type="character"),
  make_option(c("--output"), type="character"),
  make_option(c("--normalisation"), type="character",
              default= "none", help = "none/quantiles"),
  make_option(c("--extract_neighborhood"), action="store_true", default=FALSE,
              help="Extract sequence features? Takes a while..."),
  make_option(c("--smallest_unique_groups"), action="store_true", default=FALSE,
              help="Simplify protein groups"),
  make_option(c("--organisms"), type="character", default="",
              # default = "ecoli,homo_sapiens",
              help="Organisms present in the data.
              A fasta file with the same name under the fasta dir is required.
              Each organism must be provided separated by ,.
              If more than one organism is passed, the annotation_df argument is mandatory"),
  make_option(c("--filetype"), type="character", help="MaxQuant/moFF")
  
)

cat("Parsing input parameters and reading data")

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
exception <- F
# Check input arguments
if(is.null(opt$exp)) {
  cat("Please supply the experimental design with --exp")
  exception <- T
} else if((length(strsplit(opt$organisms, split = ",") %>% unlist) > 1) & is.null(opt$annotation_df)) {
  cat("Supply annotation df")
  exception <- T
} else if(!any(opt$filetype %in% c("MaxQuant", "moFF"))) {
  cat("Please supply a valid filetype: either MaxQuant or moFF")
  exception <- T
} else  if(!file.exists(opt$pepf)) {
  cat("Peptide file not found")
  exception <- T
} else if (!file.exists(opt$exp)) {
  cat("Experimental design file not found")
  exception <- T
} else if(opt$protein_file != "" & !file.exists(opt$protein_file)) {
  cat("Supplied protein file not found")
  exception <- T
}

if (exception) {
  cat("Usage:")
  quit(status=1)
}

library(GenomicRanges)
library(readxl)
library(Biostrings)
library(MSqRob)
library(MSnbase)
library(tidyr)
library(seqinr)


extract_neighborhood <- opt$extract_neighborhood
organisms <- strsplit(opt$organisms, split=",") %>% unlist
exp_annotation <- read.table(file = opt$exp, header = T, sep = "\t")
file_proteinGroups <- opt$protein_file
file_peptides_txt <- opt$pepf

if(!is.null(opt$annotation_df)) {
  annotation_df <- read.table(file = opt$annotation_df, sep = "\t", header=T, stringsAsFactors = F)
}

# exp_annotation <- read.table(file = "thp1/exp_annotation.tsv", header = T, sep = "\t")
# file_peptides_txt <- "thp1/peptide_summary_intensity_moFF_run_pellet.tab"
# opt <- list(normalisation = "quantiles", smalles_unique_grops=T, filetype="moFF")


# Import MSnSet
peptides_set <- import2MSnSet(file_peptides_txt, filetype=opt$filetype,remove_pattern = T)



###################################

cat("Preprocessing MSnSet")
# Standard preprocessing of MaxQuant data
if(opt$filetype == "MaxQuant") {
  peptides_set_preprocessed <- preprocess_MaxQuant(peptides_set, # the MSnset containing the peptides results
                                                   accession="Proteins", # column containing the protein identifiers, required if smallestUniqueGroups is TRUE
                                                   exp_annotation=exp_annotation, # experimental annotation to take batch effect into account
                                                   logtransform=TRUE, # log transform intensities (to make them more symmetrical  https://www.sciencedirect.com/science/article/pii/S1874391917301239?via%3Dihub)
                                                   base=2, # base to use in the log transform
                                                   # normalisation="quantiles", # remove systematic bias across samples (runs)
                                                   normalisation=opt$normalisation,
                                                   smallestUniqueGroups=opt$smallest_unique_groups, # occam razor for the protein groups
                                                   useful_properties=c("Proteins","Sequence"), # columns to keep in further analysis
                                                   filter=c("Potential.contaminant","Reverse"), # remove entries positive for these fields
                                                   remove_only_site=TRUE, # true if we want to discard proteins identified with modified peptides
                                                   file_proteinGroups=file_proteinGroups, # this needs to be passed if remove_only_site is True
                                                   filter_symbol="+", # if +, the field is positive
                                                   minIdentified=2 # the peptide needs to have been identified at least twice in the dataset
                                                   )
  cat("Building protdata object")
  system.time(proteins_set <- MSnSet2protdata(peptides_set_preprocessed, accession="Proteins"))
  
} else if (opt$filetype == "moFF") {
  peptides_set_preprocessed <- preprocess_MSnSet(MSnSet = peptides_set,
                                                 accession = "prot",
                                                 exp_annotation = exp_annotation,
                                                 logtransform = TRUE,
                                                 base = 2,
                                                 normalisation = opt$normalisation,
                                                 smallestUniqueGroups = opt$smallest_unique_groups,
                                                 split = ", ",
                                                 useful_properties = c("prot", "peptide"),
                                                 minIdentified = 2)
  cat("Building protdata object")
  system.time(proteins_set <- MSnSet2protdata(peptides_set_preprocessed, accession="prot"))
  
}
# Introduce into the MSqRob pipeline
# MSnSet -> protdata
cat("Generating DataFrame")
data <- getData(proteins_set)
data_df <- data %>% do.call(rbind, .)
accessions <- getAccessions(proteins_set)
data_df$Proteins <- rep(accessions, lapply(data, nrow) %>% unlist)
peptides <- data_df


## SEQUENCE FEATURE EXTRACTION
#############################################################################################

# if(opt$extract_neighborhood) {
#   cat("Extrating sequence features")
#   peptides <- peptides %>% select(Sequence, Proteins)
#   peptides$Sequence <- as.character(peptides$Sequence)
#   peptides$Proteins <- as.character(peptides$Proteins)
#   
# 
#     
#   # Read fasta files provided in organisms
#   fasta_files_list <- lapply(organisms, function(x) paste0("fasta/", x, ".fasta"))
#   fastas_list <- lapply(fasta_files_list, function(f) read.fasta(file = f, seqtype = "AA"))
#   names(fastas_list) <- organisms
#   
#   # Select just the Sequence and Proteins field
#   # Keep the first instance of each peptide
#   # There's one per run and we want to have just one
#   peptides <- peptides[(!peptides %>% duplicated),]
#   
#  
#   if(!is.null(opt$annotation_df)) {
#       peptides <- left_join(peptides, annotation_df %>% filter(!is.na(taxon)), by = "Proteins")
#   }
#   
#   # Refine the list of proteomes (fastas_list) so that the keys
#   # contain just the protein name
#   # sp|PXXXX|YYYY -> PXXXX
#   proteins_list <- lapply(fastas_list, function(fasta_file) {
#     getAnnot(fasta_file) %>%
#       lapply(., function(x) strsplit(x, split = "\\|") %>% unlist %>% .[2])
#   })
#   names(proteins_list) <- names(fastas_list)
# 
#   
#   # peptides %>% group_by(taxon) %>% summarise(count=n())
#   # peptides %>% group_by(taxon) %>% summarise(count=length(unique(Proteins)))
#   # peptides <- peptides %>% arrange(taxon, Proteins)
#   
#   # Initialize a list to store the matches and the not found 
#   match_list <- list()
#   not_found_list <- list()
#   
#   n_prots <- 1
#   k <- 1
#   
#   # Create a new prot field that will store the name of the first protein in the protein group
#   # without isoform annotation
#   peptides$prot <- peptides$Proteins %>% strsplit(., split=";") %>% lapply(., function(x) x[1]) %>%
#     lapply(., function(x) strsplit(x, split = "-") %>% unlist %>% .[1]) %>%
#     unlist
# 
#   organisms_full <- c("Escherichia coli (strain K12)", "Homo sapiens")
#   names(organisms_full) <- organisms
#   
#   # This block of code goes through the proteomes
#   # to find the sequence of each protein
#   # Takes a couple of minutes
#   ################################################################################
#   for (org in organisms) {
#     
#     if(is.null(annotation_df)) {
#       taxon_peptides <- peptides
#     } else {
#       taxon_peptides <- peptides %>% filter(taxon == organisms_full[org])
#     }
#     taxon_proteins <- proteins_list[[org]]
#     taxon_fasta <- fastas_list[[org]]
#     proteins <- taxon_peptides$prot %>% unique
#     
#     cat(paste0("Analyzing proteome ", organisms_full[org]))
#     pb <- txtProgressBar(min = 0, max = length(proteins), initial = 0) 
#     for (p in proteins) {
#       p_peptides <- taxon_peptides %>% filter(prot == p)
#       # split the protein group in the constituent protein ids and remove isoforms
#       # ps <- strsplit(protein_group, split=";") %>% unlist %>% lapply(., function(x) strsplit(x, split = "-") %>% unlist %>% .[1]) %>% unlist %>% unique
#       fasta_pos <- which(taxon_proteins %in% p)
#       if(length(fasta_pos) > 0) {
#         protein_sequence <- getSequence(taxon_fasta[fasta_pos])[[1]] %>%
#           paste(., collapse="") %>%
#           # lapply(., function(x) paste(x, collapse="")) %>% unlist
#           AAString
#         # Biostrings::(x = protein_sequences)
#         # for (pep in p_peptides$Sequence) {
#         # print(pep)
#         match_range_list <- lapply(p_peptides$Sequence, function(pep) matchPattern(pattern = pep, subject = protein_sequence))
#         match_list[[p]] <- match_range_list
#         
#         # iranges_list[[i]] <- match_range
#         # start(match_range) <- max(1, start(match_range) - 15)
#         # end(match_range) <- min(width(protein_sequences[1]), end(match_range) + 15)
#         # window_15[[pep]] <- protein_sequences[[1]][match_range]
#         # i <- i + 1
#         # } # end of peptide
#       } else {
#         # print(paste0(p, " not found"))
#         not_found_list[[k]] <- p
#         match_list[[p]] <- NULL
#         
#         k <- k + 1
#       }
#       # if (n_prots %% 500 == 0) print(paste0(n_prots, " analyzed"))
#       n_prots <- n_prots + 1
#       setTxtProgressBar(pb,n_prots)
#       
#     } # end of protein
#   } # end of taxon
#   ################################################################################
#   # match_list %>% length
#   # not_found_list %>% length
#   
#   
#   n_rows <- match_list %>% lapply(length) %>% unlist %>% sum
#   result <- matrix(nrow=n_rows, ncol=3)
#   colnames(result) <- c("pep", "neigh", "prot")
#   window_length <- 15
#   i <- 1
#   k <- 1
#   pb <- txtProgressBar(min = 0, max = length(match_list), initial = 0) 
#   # This block of code goes through the protein sequences
#   # to find the 15 aminoacid window
#   # For now, it uses the window of the first sequence (ignores all the other matches)
#   # Takes around 5 mins
#   cat("Extracting 15 aminoacid window")
#   ################################################################################
#   for(p in names(match_list)) {
#     for(m in match_list[[p]]) {
#       if (length(m) == 1) {
#         pep <- as.character(m)
#         start(m) <- max(1, start(m) - window_length);
#         end(m) <- min(m@subject %>% length , end(m) + window_length)
#         result[i, 1] <- pep
#         result[i, 2] <- as.character(m)
#         result[i, 3] <- p
#         i <- i + 1
#       } else {
#         print(paste0("Repetitions found for ", as.character(m)[1]))
#       }
#     }
#     k <- k + 1
#     setTxtProgressBar(pb,k)
#     # if(k %% 100==0) print(k)
#   }
#   
#   result_df <- result %>% as.data.frame()
#   head(result_df)
# 
# ## END
# #############################################################################################
# # Keep the peptides for which features could be extracted
# data_df_subset <- data_df %>% filter(Sequence %in% result[, "pep"])
# } else {
  data_df_subset <- data_df
# }

# Get organism link
if(!is.null(opt$annotation_df)) {
  cat("Introducing annotation")
  data_df_subset <- left_join(data_df_subset, annotation_df, by = "Proteins")
}
# Make data wide
data_df_subset_wide <- data_df_subset %>% select(-condition) %>% spread(run, quant_value)


# Remove missing data
cat("Removing missing data")
data_df_subset_wide <- data_df_subset_wide[complete.cases(data_df_subset_wide),]

if(is.null(opt$annotation_df)) {
  data_df_subset_wide <- cbind(data_df_subset_wide[,1:2], "taxon"=NA, data_df_subset_wide[,3:ncol(data_df_subset_wide)])
}


colnames(data_df_subset_wide)[1:3] <- c("peptide", "protein", "taxon")
# data_df_subset_wide %>% group_by(taxon) %>% summarise(count = n())

## Keep proteins with at least 2 peptides available
cat("Filtering peptides: minimum 2 peptides / protein")
solid_proteins <- which((data_df_subset_wide$protein %>% table) >= 2) %>% names
data_df_subset_wide_solid <- data_df_subset_wide %>% filter(protein %in% solid_proteins)

# Sort both tables in the same order!!
cat("Arranging data by peptide sequence")
data_df_subset_wide_solid$peptide <- as.character(data_df_subset_wide_solid$peptide)
data_df_subset_wide_solid <- data_df_subset_wide_solid %>% arrange(peptide)

cat(paste0("Exporting data for ", nrow(data_df_subset_wide_solid), " peptides"))

if(opt$extract_neighborhood) {
  
  result_df <- result_df[-(result_df %>% is.na %>% which(arr.ind = T) %>% .[, 1]),]
  result_df <- result_df[result_df[,"pep"] %in% data_df_subset_wide_solid$peptide,]
  result_df <- result_df %>% arrange(pep)
  all((data_df_subset_wide_solid$peptide == as.character(result_df$pep)))
  write.table(file = file.path(opt$output, "peptides_neighborhood.tsv"), sep = "\t", x = result_df, col.names=T, row.names=F, quote=F)
}

write.table(file = file.path(opt$output, "ms1_intensities.tsv"), sep = "\t",
            x = data_df_subset_wide_solid %>% select(-peptide), col.names=T, row.names=F, quote=F)

write.table(file = file.path(opt$output, "peptides.tsv"), sep = "\t",
            x = data_df_subset_wide_solid %>% select(peptide), col.names=T, row.names=F, quote=F)
