library("dplyr")
library("optparse")
library("MSnbase")

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

##moff data into msnset
## read moff file in msnset object
path <- file.path(output_moff, "peptide_summary_intensity_moFF_run_xic_norm_norm.tab")
set <- readMSnSet2(path,ecol = -c(1,2), sep = '\t')
peptide_summary_intensity <- read.table(file = path, sep = '\t', header=T, stringsAsFactors = F)

peptide_summary_intensity$Protein.IDs %>% strsplit(., split = ", ") %>% unlist

# set2 <- readMSnSet2(file.path(output_moff, "peptide_summary_intensity_moFF_run.tab"),ecol = -c(1,2), sep = '\t')
### optional

experimental_design <- read.table(file.path(root_dir, exp_name, "data", "experimental_design.tsv"), sep = "\t", header=T, stringsAsFactors = F) %>% filter(Fraction==1)
## meta data for each sample; dataframe with one row/sample
# pd = data.frame(condition = experimental_design$Experiment, replicate =  experimental_design$Replicate)

# this assumes that the sampleNames in set i.e the column names of the peptide_summary file
# are the same, though not necessarily in the same order
# as the result of pasting the fields Experiment and Replicate from the experimental_design.tsv

# rownames(pd) = paste(pd$condition, pd$replicate, sep = "_")
# pd <- pd[sampleNames(set),]

pd <- do.call(rbind, colnames(peptide_summary_intensity[, -(1:2)]) %>% strsplit(., split = "_")) %>% as.data.frame
colnames(pd) <- c("condition", "replicate")
rownames(pd) <- colnames(peptide_summary_intensity[, -(1:2)])


## meta data for each peptide; dataframe with one row/peptide
fd = peptide_summary_intensity %>% select(Protein.IDs, Sequence)
rownames(fd) = featureNames(set)
## add meta data to set object
set = MSnSet(exprs(set),
             fData =  AnnotatedDataFrame(fd),
             pData = AnnotatedDataFrame(pd))

##robust summarisation
protset <- combineFeatures(set, fun="median", groupBy = fData(set)$Protein.IDs, cv = FALSE)

quantification <- exprs(protset)
write.table(x = quantification,
            file = file.path(root_dir, exp_name, "quantification", "msnbase_quantification.tsv"),
            col.names = T, row.names = T, quote = F, sep = "\t")

# protein_ids <- rownames(quantification) %>% strsplit(., split = ", ") %>% unlist
# protein_ids[protein_ids %>% duplicated]
# table(protein_ids) %>% sort %>% rev %>% head(100)
# grep("Q9NQ29", rownames(quantification), value =T)

