library("dplyr")
library("ggplot2")
library("gProfileR")
library("optparse")

home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
option_list <- list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "thesis", "genedata")),
  make_option(c("--exp_name"), type="character", default="thp1"),
  make_option(c("--suffix"), type="character", default=""),
  make_option(c("--log2fc_threshold"), default=1)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
suffix <- opt$suffix
log2fc_threshold <- opt$log2fc_threshold

if(exp_name == "maxlfq") source(file.path(root_dir, "scripts", "R", "check_organism.R"))

input_dir <- file.path(root_dir, exp_name, "quantification")
RSqM_signif_file <- file.path(input_dir, paste0("RSqM_signif", suffix, ".tsv"))
RSqM_signif <- read.table(file = RSqM_signif_file, header = T, sep = "\t", stringsAsFactors = F)
RSqM_signif <- RSqM_signif[-(RSqM_signif$Protein.IDs %>% grep(pattern = "CONTAMINANT", x = .)),]

if(exp_name == "maxlfq") {
  organism_check <- check_organism(as.character(RSqM_signif$Protein.IDs), split = ", ")
  RSqM_signif$taxon <- organism_check[[2]]
  RSqM_signif <- RSqM_signif %>% filter(!is.na(taxon))
  write.table(x = RSqM_signif, file = file.path(root_dir, exp_name, "quantification", "RSqM_signif_annotated.tsv"), quote = F,row.names = F,col.names = T,sep = "\t")
}


plots_dir <- file.path(root_dir, exp_name, "plots")

volcano_plot <- ggplot(data = RSqM_signif, aes(x = estimate, y = -log10(qval))) +
  geom_vline(xintercept = c(-1, 1)) +
  geom_hline(yintercept = -log10(0.05))

if(exp_name == "maxlfq") {
  volcano_plot <- volcano_plot + geom_point(aes(col = taxon))
} else {
  volcano_plot <- volcano_plot + geom_point()
}

volcano_plot

ggsave(filename = file.path(plots_dir, paste0("volcano_plot", suffix, ".png")), plot = volcano_plot)

histogram <- ggplot(data = RSqM_signif, aes(x = estimate)) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed")

if(exp_name == "maxlfq") {
  histogram <- histogram + geom_histogram(bins = 100, position="dodge", aes(fill=taxon))
} else {
  histogram <- histogram + geom_histogram(bins = 100, position="dodge")
}

histogram


ggsave(filename = file.path(plots_dir, paste0("histogram", suffix, ".png")), plot = histogram)

if(exp_name != "maxlfq") {
 
 up_prots <- RSqM_signif %>% filter(signif) %>% filter(estimate > log2fc_threshold) %>% .$Protein.IDs %>% strsplit(split = ", ") %>% unlist
 down_prots <- RSqM_signif %>% filter(signif) %>% filter(estimate < -log2fc_threshold) %>% .$Protein.IDs %>% strsplit(split = ", ") %>% unlist
 stable_prots <- RSqM_signif %>% filter(signif) %>% filter(abs(estimate) < log2fc_threshold) %>% .$Protein.IDs %>% strsplit(split = ", ") %>% unlist
 dap_prots <- c(up_prots, down_prots)
 na_prots <- RSqM_signif %>% filter(is.na(estimate)) %>% .$Protein.IDs
 
 write.table(x = up_prots, file = file.path(input_dir, "up_regulated.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = down_prots, file = file.path(input_dir, "down_regulated.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = dap_prots, file = file.path(input_dir, "dap_prots.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = stable_prots, file = file.path(input_dir, "stable_prots.txt"), quote = F, row.names = F, col.names = F)
 write.table(x = na_prots, file = file.path(input_dir, "na_prots.txt"), quote = F, row.names = F, col.names = F)

 # RSqM_signif %>% filter(estimate == min(estimate, na.rm = T))
 
 gsea_estimated <- gprofiler(query = c(dap_prots, stable_prots), organism = "hsapiens",
                             png_fn = file.path(input_dir, "gsea.png"),
                             include_graph = T)
 
 gsea_not_estimated <- gprofiler(query = c(na_prots), organism = "hsapiens")
 
 write.table(x = gsea_estimated, file = file.path(input_dir, "gsea_estimated.tsv"), sep = "\t", quote = F, row.names = F)
 write.table(x = gsea_estimated %>% select(term.id, p.value), file = file.path(input_dir, "gsea_estimated_term.txt"),
             sep = "\t", quote = F, row.names = F, col.names = F)
 
 
 library("UniProt.ws")
 library("org.Hs.eg.db")
 library("annotate")
 up <- readRDS(file = file.path(root_dir, exp_name, "quantification", "uniprot_ws.rds"))
 keys <- stable_prots
 columns <- c("UNIPROTKB", "ENTREZ_GENE")
 kt <- "UNIPROTKB"
 res <- UniProt.ws::select(x = up, keys = keys, columns = columns, keytype = kt)
 entrez_gene_names <- getSYMBOL(res$ENTREZ_GENE %>% na.omit, data ='org.Hs.eg')
 res <- left_join(data.frame(ENTREZ_GENE = na.omit(res$ENTREZ_GENE), name = entrez_gene_names), res, by = "ENTREZ_GENE")

 webapp <- read.table(file.path(root_dir, exp_name, "quantification", "uniprot_id_mapping.txt"),
                      header=T, col.names = c("UNIPROTKB", "webapp_name"))
 
 res <- full_join(res, webapp, by = "UNIPROTKB") %>% select(UNIPROTKB, name, webapp_name)
 # Generate multifasta file for KEGG KAAS
 # library("seqinr")
 # fasta_database <- read.fasta(file = file.path(root_dir, exp_name, "databases", "all.fasta"), seqtype = "AA")
 # annotations <- seqinr::getAnnot(fasta_database)
 # 
 # indices <- lapply(dap_prots, function(x) grep(pattern = x, x = annotations)) %>% unlist
 # fasta_subset <- fasta_database[indices]
 # write.fasta(sequences = getSequence(fasta_subset), names = getAnnot(fasta_subset),
 #             file.out = file.path(root_dir, exp_name, "kegg", "kegg_multifasta.fasta"))
 write.table(x = res$name, file = file.path(root_dir, exp_name, "quantification", "compath.txt"), row.names = F, col.names = F, quote = F)
 # http://compath.scai.fraunhofer.de/query
}
