library("optparse")
option_list <- list(
  make_option(c("--peptide_file"), type="character",
              default="MaxLFQ/peptides.txt"),
  make_option(c("--output_dir"), type="character",
              default=".")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
peptide_file <- opt$peptide_file

source("check_organism.R")
proteins <- unique(read.table(file = peptide_file, header=T,sep="\t",stringsAsFactors = F)$Proteins)
organism_check <- check_organism(proteins = proteins, split = ";")
annotation_df <- data.frame(Proteins = proteins, taxon = organism_check[[2]]) %>% na.omit
write.table(x = annotation_df,
            file = file.path(opt$output_dir, "annotation_df.tsv"), sep="\t",
            quote = F,row.names = F,col.names = T
            )
