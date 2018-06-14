#' Smallest unique protein groups
#'
#' @description For a given vector of protein group names, outputs the names of those protein groups for which none of its member proteins is present in a smaller protein group.
#' @param proteins A vector of characters or factors containing single proteins and/or protein groups (i.e. proteins separated by a separator symbol).
#' @param split The character string that is used to separate the indivudual protein names in each protein group.
#' @return A character vector containing the names of the protein groups for which none of its proteins is present in a smaller protein group.
#' @examples #This example will give the names of the protein groups in the MSnSet object peptidesCPTAC for which none of its proteins is present in a smaller protein group.
#' library("MSnbase")
#' #Select the columns containing intensities
#' colInt <- grepEcols(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), pattern="Intensity.", split = "\t")
#' #Import MaxQuant's peptides.txt file and convert it to an MSnSet object
#' peptidesCPTAC <- readMSnSet2(system.file("extdata/CPTAC", "peptides.txt", package = "MSqRob"), ecol = colInt, sep = "\t")
#'
#' #Our approach: a peptide can map to multiple proteins, as long as there is none of these proteins present in a smaller subgroup.
#' groups <- smallestUniqueGroups(fData(peptidesCPTAC)$Proteins)
#' groups
#' @export
#' 
# smallestUniqueGroups <- function(proteins, split=";"){
# 
#   # list of unique protein groups where every element is represented by a character vector
#   # where every element is a character string representing a protein id
#   b <- strsplit(x=as.character(unique(proteins)), split=split, fixed=TRUE)
#   # original_b <- b
# 
#   erbij <- vector()
# 
#   j <- 1
#   while(length(b)!=0)
#   {
#     # erbij is a character vector storing protein groups
#     # each protein group is represented as a string (not a character vector) with each protein id being separated by "split"
#     # when j = n, erbij contains protein groups of length < n that have pass the Occam razor filter
#     # After this line of code, it gets added those of length n
#     # if j = 1, prior to running this code, erbij is empty
#     
#     erbij <- c(erbij, sapply(b[sapply(b, length)==j], function(x) paste(x, collapse=split)))
# 
#     # a contains the protein ids available in groups of length j in a character string 
#     a <- unlist(b[sapply(b, length) == j])
#     # b is a list with the protein groups of length greater than j
#     # these are the protein groups that we can apply the Occam's razor on.
#     b <- b[sapply(b, length)>j]
# 
#     if(length(b) != 0){
#     # welke will store T/F for every protein group f length > j
#     welke <- vector()
#     for(i in 1:length(b))
#     {
#       # welke_i is true if none of the ids in b_i is was seen in protein groups of length == j
#       # welke_i is false if at least one of the ids in b_i is was seen in protein groups of length == j
#       welke[i] <- !any(b[[i]] %in% a)
#     }
#     # select the protein groups that contain unseen ids
#     b <- b[welke]
#     j <- j+1
#     }
#   }
# 
#   erbij <- unlist(erbij)
#   return(erbij)
# }


smallestUniqueGroups <- function(proteins, split=";"){
  
  b <- strsplit(x=as.character(unique(proteins)),split=split,fixed=TRUE)
  
  erbij <- vector()
  
  j <- 1
  while(length(b)!=0)
  {
    
    erbij <- c(erbij,sapply(b[sapply(b, length)==j], function(x) paste(x, collapse=split)))
    
    a <- unlist(b[sapply(b, length)==j])
    b <- b[sapply(b, length)>j]
    
    if(length(b)!=0){
      welke <- vector()
      for(i in 1:length(b))
      {
        welke[i] <- !any(b[[i]] %in% a)
      }
      
      b <- b[welke]
      j <- j+1
    }
  }
  
  erbij <- unlist(erbij)
  return(erbij)
}

library("optparse")
library("dplyr")

option_list = list(
  make_option(c("--output_moff"), type="character", default="//hest/aoj//thesis/genedata/maxlfq/peptideShaker_out/PSM_reports/output_moff_RAW/")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
output_moff <- opt$output_moff

peptide_summary_intensity <- read.table(file = file.path(output_moff, "peptide_summary_intensity_moFF_run.tab"), header=T, sep = "\t")
proteins <- peptide_summary_intensity$prot
smallest_groups <- smallestUniqueGroups(proteins, split = ", ")

###################################################

# smallest_groups %>% sapply(length) %>% table
# smallest_groups %>% duplicated %>% table
# smallest_groups

###################################################

peptide_summary_intensity_occam <- filter(peptide_summary_intensity, prot %in% smallest_groups)

x <- peptide_summary_intensity_occam$prot %>% as.character %>% unique %>% strsplit(., split = ", ") %>% unlist %>% table
x[x != 1]

write.table(x = peptide_summary_intensity_occam,
            file = file.path(output_moff, "peptide_summary_intensity_moFF_run_occam.tab"),
            col.names = T, row.names = F, quote = F, sep = "\t")
