library("dplyr")
library("ggplot2")
library("gProfileR")
library("optparse")
library(kableExtra)
library(xtable)
source("read_quant.R")

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "thesis", "genedata")),
  make_option(c("--exp_name"), type="character", default="thp1"),
  make_option(c("--input_dir"), type="character", default="quantification"),
  make_option(c("--output_dir"), type="character", default="."),
  make_option(c("--suffix"), type="character", default=""),
  make_option(c("--quant_file"), type="character", default=""),
  make_option(c("--log2fc_threshold"), default=1)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
suffix <- opt$suffix
plots_dir <- file.path(opt$output_dir, "plots")
log2fct <- opt$log2fc_threshold

input_dir <- opt$input_dir

if(interactive()) opt <- list(log2fc_threshold = 1, quant_file = "../thp1/RSqM_signif.tsv",
                              output_dir = "/home/antortjim/MEGA/Master/Thesis/Report", suffix="_msqrob")


quant <- read_quant(filename = opt$quant_file, filetype = "MSqRob")
RSqM_signif <- quant
RSqM_signif %>% filter(signif) %>% nrow()

nrow(quant)

volcano_plot <- ggplot(data = quant, aes(x = estimate, y = -log10(qval))) +
  geom_vline(xintercept = c(-log2fct, log2fct), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_point() +
  labs(x="log2FC")



ggsave(filename = file.path(plots_dir, paste0("volcano_plot", opt$suffix, ".png")), plot = volcano_plot)
ggsave(filename = file.path("../../Report/plots/", paste0("volcano_plot", opt$suffix, ".png")), plot = volcano_plot)

histogram <- ggplot(data = quant, aes(x = estimate)) +
  geom_vline(xintercept = c(-log2fct, log2fct), linetype="dashed") +
  geom_histogram(bins = 100, position="dodge", fill = "lightgreen", col = "#aaaaaa") +
  labs(x="log2FC", y = "Count")


ggsave(filename = file.path(plots_dir, paste0("histogram", opt$suffix, ".png")), plot = histogram)
ggsave(filename = file.path("../../Report/plots/", paste0("histogram", opt$suffix, ".png")), plot = histogram)

 up_prots <- quant %>% filter(signif) %>% filter(estimate > log2fct) %>% .$Protein.IDs
 down_prots <- quant %>% filter(signif) %>% filter(estimate < -log2fct) %>% .$Protein.IDs
 stable_prots <- quant %>% filter(signif) %>% filter(abs(estimate) < log2fct) %>% .$Protein.IDs
 signif_prots_group <- list(up_prots, down_prots, stable_prots)
 length(unlist(signif_prots_group))
 
 signif_prots_list <- signif_prots_group %>% lapply(function(x) strsplit(x, split = ", ") %>% unlist)
 up_prots <- signif_prots_list[[1]]
 down_prots <- signif_prots_list[[2]]
 stable_prots <- signif_prots_list[[3]]
 dap_prots <- c(up_prots, down_prots)
 na_prots <- quant %>% filter(is.na(estimate)) %>% .$Protein.IDs
 signif_prots <- c(dap_prots, stable_prots)
 write.table(x = up_prots, file = file.path(opt$output_dir, "up_regulated.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = down_prots, file = file.path(opt$output_dir, "down_regulated.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = dap_prots, file = file.path(opt$output_dir, "dap_prots.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = stable_prots, file = file.path(opt$output_dir, "stable_prots.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = na_prots, file = file.path(opt$output_dir, "na_prots.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = signif_prots_group %>% unlist, file = file.path(opt$output_dir, "signif_prots.txt"), quote = F, row.names = F, col.names = F)
 
 # row_ids <- lapply(signif_prots, function(x) grep(pattern = x, x = quant$Protein.IDs)) %>% unlist %>% unique
  # report_table <- quant[row_ids,] %>% .[,c("Protein.IDs", "estimate", "qval")] %>% na.omit
 signif_prots <- grep(pattern = ", ", x = unlist(signif_prots_group), invert = T, value=T)
 row_ids <- which(quant$Protein.IDs %in% signif_prots)
 
 report_table <- quant[row_ids, c("Protein.IDs", "estimate", "qval")] %>%
   na.omit
 colnames(report_table)[1:2] <- c("Protein", "log2FC")
 # report_table <- cbind(Protein = report_table[,1], round(report_table[,2:3], digits = 3))
 library("UniProt.ws")
 library("org.Hs.eg.db")
 library("annotate")
 up <- UniProt.ws(taxId=9606)
 keys <- signif_prots
 columns <- c("UNIPROTKB", "ENTREZ_GENE", "PROTEIN-NAMES","ENTRY-NAME")
 kt <- "UNIPROTKB"
 UniProt_results <- UniProt.ws::select(x = up, keys = keys, columns = columns, keytype = kt)
 entrez_gene_names <- getSYMBOL(UniProt_results$ENTREZ_GENE %>% na.omit, data ='org.Hs.eg')
 res <- left_join(data.frame(ENTREZ_GENE = na.omit(UniProt_results$ENTREZ_GENE), name = entrez_gene_names),
                  UniProt_results, by = "ENTREZ_GENE")

res <- res[,c("UNIPROTKB", "PROTEIN-NAMES")]
colnames(res) <- c("Protein", "Name") 

nrow(report_table)
report_table <- left_join(report_table, res, by = "Protein") %>% dplyr::select(
  Protein, Name, log2FC, qval
)
 
report_table$Name <- report_table$Name %>% strsplit(., split = " \\(") %>% lapply(function(x) x[1]) %>% unlist

rws <- seq(1, (nrow(report_table)-1), by = 2)
col <- rep("\\rowcolor[gray]{0.95}", length(rws))

xtable(report_table,
       digits=-2,
       booktabs = TRUE,
       add.to.row = list(pos = as.list(rws)),
       command=col)
 

gsea_estimated <- gprofiler(query = grep(pattern = ", ", invert = T, value = T, x = unlist(signif_prots_group)),
          organism = "hsapiens"
          # png_fn = file.path(plots_dir, "gsea.png"),
          # include_graph = T
          )

report_table <- gsea_estimated %>% arrange(p.value) %>%
  dplyr::select(term.id, term.name, query.size, term.size, overlap.size, p.value)
colnames(report_table) <- c("Term ID", "Term Name", "Query", "Term", "Ov", "P-value")
report_table$`P-value` <- report_table$`P-value` %>% formatC(., format = "e", digits = 2)
report_table$`Term Name` <- report_table$`Term Name` %>% strsplit(., split = "\\[") %>% lapply(., function(x) x[1]) %>% unlist
rws <- seq(1, (nrow(report_table)-1), by = 2)
col <- rep("\\rowcolor[gray]{0.95}", length(rws))

print.xtable(xtable(report_table %>% head(20),
                    booktabs = TRUE,
                    add.to.row = list(pos = as.list(rws)),
                    command=col),
             type="latex", floating = F, file = file.path(opt$output_dir, "tables/gsea.tex"))



gsea_estimated$intersection[c(5,7,14,15,19)] %>% strsplit(., split=",") %>% unlist %>% unique


 write.table(x = gsea_estimated, file = file.path(opt$output_dir, "gsea_estimated.tsv"), sep = "\t", quote = F, row.names = F)
 write.table(x = gsea_estimated %>% select(term.id, p.value), file = file.path(opt$output_dir, "gsea_estimated_term.txt"),
             sep = "\t", quote = F, row.names = F, col.names = F)
 
 
# kegg <- read.table(file = "../thp1/KEGG Enriched Pathways.csv", sep = ",", header=T) %>% .[,c(2,4,5,6)]
# reactome <- read.table(file = "../thp1/Reactome Enriched Pathways.csv", sep = ",", header=T) %>% .[,c(2,4,5,6)]
# pathways <- rbind(cbind(head(kegg, 10), database="kegg"), cbind(head(reactome, 10), database="reactome"))
# pathways$Pathway.Name <- pathways$Pathway.Name %>% as.character %>% strsplit(., split = " \\- ") %>% lapply(function(x) x[1]) %>% unlist
# colnames(pathways) <- c("Pathway", "q-val", "Ov", "Size", "DB")
# pathways <- pathways %>% arrange(`q-val`)
# pathways$`q-val` <- pathways$`q-val` %>% formatC(., format = "e", digits = 2)
# xtable(pathways)
