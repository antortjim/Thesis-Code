library("dplyr")
library("optparse")
library("ggplot2")

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character"),
  make_option(c("--pipeline"), type="logical", default = F)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
pipeline <- opt$pipeline

input_dir <- file.path(root_dir, exp_name, ifelse("input_dir" %in% names(opt), opt$input_dir, "quantification"))
output_dir <- file.path(root_dir, exp_name, ifelse("output_dir" %in% names(opt), opt$output, "quantification"))



######################################################################
## Read data
######################################################################
experimental_design <- read.table(file = file.path(root_dir, exp_name, "data", "experimental_design.tsv"), header=T)
fractions_available <- F
if(any(experimental_design$Fraction != 1)) fractions_available <- T


## A from MaxQuant paper published results
###############################################################

if(!pipeline) {
  
  # The intensities file from MaxLFQ paper was generated using the script: /thesis/genedata/maxlfq/paper/load_maxlfq_data.R
  intensities <- read.table(file = file.path(input_dir, "maxlfq_peptide_intensities.tsv"), sep = "\t", header=T, check.names = F)
  peptides_Protein.IDs <- read.table(file = file.path(input_dir, "peptides_Protein.IDs.txt"), check.names = F, stringsAsFactors = F)[,1]
  peptides_sequence <- read.table(file = file.path(input_dir, "peptides_sequence.txt"), check.names = F, stringsAsFactors = F)[,1]
  combinations <- read.table(file = file.path(input_dir, "combinations.tsv"), sep = "\t") %>% as.matrix %>% .[-1,]
  combinations_names <- read.table(file = file.path(input_dir, "combinations_names.tsv"), sep = "\t", stringsAsFactors = F)[,1]
  Protein.IDs <- read.table(file = file.path(input_dir, "Protein.IDs.tsv"), sep = "\t", header=F, stringsAsFactors = F)[,1]
  
  # intensities <- intensities[,7:12]
  # combinations <- combinations[,(66-14):66]
  # combinations_names <- combinations_names[(66-14):66]
  
  ################################################################
} else {
  ## From our pipeline
  #########################################################
  #########################################################
  
  # Read moFF output
  xics <- read.table(file = file.path(input_dir, "../peptideShaker_out/PSM_reports/output_moff_RAW/peptide_summary_intensity_moFF_run.tab"),
                     sep = "\t", header=T, check.names = F, stringsAsFactors = F)
  colnames(xics)[1:2] <- c("Sequence", "Protein.IDs")
  peptides <- xics[, c("Sequence", "Protein.IDs")]
  colnames(xics) <- colnames(xics) %>% gsub(pattern = "sumIntensity_", replacement = "")
  Protein.IDs <- peptides$Protein.IDs %>% unique
  peptides$index <- 1:nrow(peptides)
  
  ## Post process moFF output so that protein ids are not shared across protein groups.
  peptides$Size <- peptides$Protein.IDs %>% strsplit(., split = ", ") %>% lapply(length) %>% unlist
  peptides <- peptides %>% arrange(-Size)
  peptides_unique <- peptides[peptides$Size == 1, ]
  
  ### SIMPLE OCCAMS RAZOR: REMOVE PEPTIDES WHERE THE PROTEIN GROUP CONTAINS A PROTEIN AVAILABLE IN A SMALLER GROUP
  ### AS STATED IN:
  ### Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics
  ### It's too slow!!!! For now we just keep protein groups of size 1
  
  # peptides_non_unique <- peptides[peptides$Size > 1, ] 
  # 
  # protein_groups_all <- peptides$Protein.IDs %>% unique %>% strsplit(., split = ", ")
  # protein_groups_unique <- peptides_unique$Protein.IDs %>% unique %>% strsplit(., split = ", ")
  # protein_groups_non_unique <- peptides_non_unique$Protein.IDs %>% unique %>% strsplit(., split = ", ")
  # 
  # group_is_subset <- function(group1, group2) {
  #   return(lapply(group1, function(x) x %in% group2) %>% unlist %>% any)
  # }
  # 
  # group_is_subset(c("A", "B"), c("E", "C", "D"))
  # 
  # contained <- logical(length = length(protein_groups_non_unique))
  # which_contained <- integer(length = length(protein_groups_non_unique))
  # protein_group_size <- peptides$Size
  # 
  # 
  # # For every group
  # for(j in 1:length(protein_groups_non_unique)) {
  #   print(j)
  #   # For every group
  #   for (i in 1:length(protein_groups_all)) {
  #     # If the groups being compared are not the same, and the second group
  #     # has not been flagged as containing proteins available in a smaller group
  #     ################
  #     # if(j != i) {
  #     #   result <- group_is_subset(protein_groups_all[[i]], protein_groups_non_unique[[j]])
  #     #   if(result) {
  #     #     contained[j] <- result
  #     #     which_contained[j] <- i
  #     #     break()
  #     #   }
  #     # }
  #     ################
  #     if(j != i) {
  #       result <- group_is_subset(protein_groups_all[[i]], protein_groups_non_unique[[j]])
  #       if(result) {
  #         contained[j] <- result
  #         break()
  #       }
  #     }
  #   }
  # }
  # 
  # # Groups that contain a protein present in a smaller group
  # which(contained)
  # # The ith element in this vector shows the id of the small group that contains a protein shared with the ith group bigger group
  # which_contained[contained]
  # 
  # peptides2 <- peptides[!contained, ]
  ################################################
  ## Keep protein groups of size 1
  ################################################
  xics <- xics[peptides$Size == 1,]
  peptides <- peptides_unique
  
  if(fractions_available) {
    write.table(x = xics, file = file.path(output_dir, "xics.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
    # call Levenberg-Marquardt optimisation
    intensities <- read.table(file = file.path(output_dir, "lm_intensities.tsv"), header=T, sep = "\t")
  }
  
  #########################################################
  #########################################################
}

# combinations <- colnames(intensities) %>% sort %>% combn(., m = 2)
combinations_names <- combinations %>% apply(., 2, function(x) paste(x, collapse = "/"))



######################################################################
## Compute peptide ratios from peptide intensity
######################################################################

peptide_ratios <- combinations %>% apply(X = ., MARGIN = 2, FUN = function(x) intensities[, x[1]] / intensities[, x[2]])
colnames(peptide_ratios) <- combinations_names
rownames(peptide_ratios) <- peptides_sequence
# if a peptide ratio is infinite, the denominator was 0
# if a peptide ratio is NaN, the denominator and the numerator were 0
# if a peptide ratio is 0, the numerator was 0
peptide_ratios[is.na(peptide_ratios) | is.infinite(peptide_ratios)] <- 0


######################################################################
## Peptides per protein
######################################################################
# peptides %>% group_by(Protein.IDs) %>% summarise(count = n()) %>% ggplot(aes(x = count)) + geom_histogram(bins = 40)
# peptides %>% group_by(Protein.IDs) %>% summarise(count = n()) %>% .$count %>% table %>% .[-1] %>% sum


######################################################################
## Compute protein ratios from peptide ratios
######################################################################

protein_ratios <- matrix(NA, nrow = length(Protein.IDs), ncol = length(combinations_names))
colnames(protein_ratios) <-  combinations_names
rownames(protein_ratios) <- Protein.IDs
status <- character(length = length(Protein.IDs))
supporting <- integer(length = length(Protein.IDs))
local_supporting <- matrix(nrow = length(Protein.IDs), ncol = ncol(combinations))
minimal_shared_peptides <- 2
i <- 1
for (protein in Protein.IDs) {
  print(i)
  # Get the peptides associated with protein protein
  protein_pos <- which(peptides_Protein.IDs == protein)
  supporting[i] <- length(protein_pos)
  # If the number of associated peptides is NOT below the threshold, 
  if(!supporting[i] < minimal_shared_peptides) {
    
    # Fetch the peptide ratios for the current protein
    protein_peptide_ratios <- peptide_ratios[protein_pos,] %>% as.matrix()
    # if("GYPHWPAR" %in% rownames(protein_peptide_ratios)) browser()
    
    # The protein ratio for comparison k is the median of the peptide ratios in comparison k 
    r <-  protein_peptide_ratios %>% apply(X = ., MARGIN = 2, FUN = function(x) median(x, na.rm = T))
    # How many peptides support each median? This number goes from 0 to the total number of peptides for this protein
    # It's the count of non-zero peptides
    number_peptides_supporting_median <- length(protein_pos) - (protein_peptide_ratios %>% apply(X = ., MARGIN = 2, FUN = function(x) x==0) %>% colSums)
    local_supporting[i,] <- number_peptides_supporting_median
    r[number_peptides_supporting_median < minimal_shared_peptides] <- 0
    status[i] <- "SUPPORT"
    # else, set the protein ratios equal to 0
  } else {
    r <- rep(0, length(combinations_names))
    status[i] <- "FRAGILE"
  }
  
  # Set the ratios of the ith protein in the ith row of the matrix
  protein_ratios[i,] <- r
  
  i <- i + 1
}

supporting %>% table


######################################################################
## Visualize distribution of the number of ratios equal to 0
## (max equal to the number of combinations, i.e 15) 
######################################################################
x <- rowSums(protein_ratios == 0)
table(x == ncol(protein_ratios))
table(x)
hist(x, breaks = 30)
# Hopefully most proteins have 0 ratios = 0, i.e the histogram should have a high bin at 0 and very low bins elsewhere
# and lower the more to the right we move in the x axis

######################################################################
## Export the protein ratios to a file
######################################################################
write.table(x = protein_ratios, file = file.path(output_dir, "protein_ratios.tsv"), row.names = T, col.names = T, sep = "\t", quote = F)
write.table(x = rownames(protein_ratios), file = file.path(output_dir, "protein_ratios.IDs.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
