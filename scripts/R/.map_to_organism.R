library("dplyr")
library("optparse")
library("seqinr")
option_list = list(
  make_option(c("--root_dir"), type="character", default="//hest/aoj//thesis/genedata"),
  make_option(c("--exp_name"), type="character", default="maxlfq")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name

output_moff <- file.path(root_dir, exp_name, "peptideShaker_out/PSM_reports/output_moff_RAW")
output_dir <- file.path(root_dir, exp_name, "quantification")

path <- file.path(output_moff, "peptide_summary_intensity_moFF_run_xic_norm_norm.tab")

peptide_summary_intensity <- read.table(file = path, sep = '\t', header=T, stringsAsFactors = F)
peptide_summary_intensity$Protein.IDs %>% strsplit(., split = ", ") %>% unlist

databases_path <- file.path(root_dir, exp_name, "databases")

databases <- list.files(path = databases_path, "uniprot")
databases_standard <- c("Escherichia coli (strain K12)", "Homo sapiens")
names(databases_standard) <- databases
id_mapping <- data.frame(NULL)
for(i in 1:length(databases)) {
  IDs <- getAnnot(read.fasta(file = file.path(databases_path, databases[i]), seqtype = "AA")) %>%
    lapply(function(x) strsplit(x, split = "\\|")[[1]][2]) %>% unlist
  
  id_mapping <- rbind(id_mapping,
                      cbind(Protein.IDs = IDs,
                            DB = databases[i])
                      )
}

id_mapping$DB <- databases_standard[id_mapping$DB]


check_id_taxonomy <- function(x, id_mapping, split = ", ") {
  ids <- strsplit(as.character(x), split = split) %>% unlist
  database <- filter(id_mapping, Protein.IDs %in% ids) %>% .$DB %>% as.character
  unambiguous <- length(unique(database)) == 1
  result <- ifelse(unambiguous, database[1], NA)
  return(result)
}
