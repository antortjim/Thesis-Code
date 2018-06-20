library(readxl)
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")

supplementary_file <- file.path(home_dir, "thesis", "genedata", "maxlfq", "mcp.M113.031591-1.xlsx")

supplementary <- read_xlsx(path = supplementary_file)
taxonomy <- supplementary[,c("Protein IDs", "Taxonomy")]
protein_ids <- taxonomy$`Protein IDs` %>% strsplit(., split = ";")
group_length <- lapply(protein_ids, length) %>% unlist

taxonomy <- data.frame(Protein.IDs = unlist(protein_ids), Taxonomy = as.character(rep(taxonomy$Taxonomy, group_length)))
check_organism <- function(proteins, split) {
  validated <- logical(length = length(proteins))
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
