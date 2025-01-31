library("readxl")
library("dplyr")
# proteins <- read.table(file = "proteins.txt", header = F, stringsAsFactors = F)[,1]
supplementary_file <- "MaxLFQ/mcp.M113.031591-1.xlsx"

supplementary <- read_xlsx(path = supplementary_file)
taxonomy <- supplementary[,c("Protein IDs", "Taxonomy")]
protein_ids <- taxonomy$`Protein IDs` %>% strsplit(., split = ";")
group_length <- lapply(protein_ids, length) %>% unlist

taxonomy <- data.frame(Protein.IDs = unlist(protein_ids), Taxonomy = as.character(rep(taxonomy$Taxonomy, group_length)))
check_organism <- function(proteins,split=";") {
  validated <- rep(F, length(proteins))
  result <- character(length = length(proteins))
  for (j in 1:length(proteins)) {
    value <- proteins[j]
    uniprot_ids <- strsplit(value, split = split) %>% unlist
    organisms <- taxonomy %>% filter(Protein.IDs %in% uniprot_ids) %>% .$Taxonomy %>% as.character
    validated[j] <- all(organisms == organisms[1])
    result[j] <- organisms[1]
    print(j)
  }
  return(list(validated, result))
}

