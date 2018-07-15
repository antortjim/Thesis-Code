library("readxl")
# proteins <- read.table(file = "proteins.txt", header = F, stringsAsFactors = F)[,1]

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj",
                   ifelse(Sys.info()["user"] == "aoj", "/z/home/aoj",
                          "/home/antortjim/"))


data_dir <- ifelse(Sys.info()["user"] == "aoj", "thesis/genedata/maxlfq/",
                   "MEGA/Master/Thesis/Code/scripts/data")

exp_dir <-  ifelse(Sys.info()["sysname"] == "Windows",
                   file.path(root_dir, exp_name),
                   file.path(home_dir, data_dir))

supplementary_file <- ifelse(Sys.info()["sysname"] == "Windows",
       "MaxLFQ/mcp.M113.031591-1.xlsx",
       file.path(exp_dir, "mcp.M113.031591-1.xlsx")
       )


supplementary <- read_xlsx(path = supplementary_file)
taxonomy <- supplementary[,c("Protein IDs", "Taxonomy")]
protein_ids <- taxonomy$`Protein IDs` %>% strsplit(., split = ";")
group_length <- lapply(protein_ids, length) %>% unlist

taxonomy <- data.frame(Protein.IDs = unlist(protein_ids), Taxonomy = as.character(rep(taxonomy$Taxonomy, group_length)))
check_organism <- function(proteins,split=";") {
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

