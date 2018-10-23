home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
library("dplyr")
library("ggplot2")
library("VennDiagram")
root_dir <- file.path("thesis", "genedata")
exp_name <- "thp1"
output_dir <- file.path(home_dir, root_dir, exp_name, "benchmark")
genedata <- read.table(file = file.path(home_dir, root_dir, exp_name, "genedata_processed", "volcano_data.tsv"), header=T, stringsAsFactors = F)
colnames(genedata)[colnames(genedata) == "name"] <- "Protein.IDs"

pipeline <- read.table(file = file.path(home_dir, root_dir, exp_name, "quantification", "RSqM_signif.tsv"), header=T, sep = "\t", stringsAsFactors = F) %>%
  filter(!is.na(estimate))

# Remove results of the pipeline for which the estimate is 0
# pipeline <- pipeline %>% filter(abs(estimate) > 0.1)
  
benchmark_dataset <- data.frame(
  Protein.IDs = c(pipeline$Protein.IDs, genedata$Protein.IDs),
  estimate = c(pipeline$estimate, genedata$log2FC),
  signif = c(pipeline$signif, genedata$criteria %in%  c("PVAL", "BOTH")),
  workflow = rep(c("pipeline", "genedata"), c(nrow(pipeline), nrow(genedata)))
  )

benchmark_dataset %>% group_by(workflow, signif) %>% summarise(count = n())

benchmark_dataset <- benchmark_dataset[-grep(", ", benchmark_dataset$Protein.IDs),]
benchmark_dataset <- benchmark_dataset[-grep("CONTAMINANT", benchmark_dataset$Protein.IDs),]

benchmark_dataset %>% group_by(workflow, signif) %>% summarise(count = n())

x <- lapply(benchmark_dataset$workflow %>% unique, function(w) {
  lapply(c(T, F), function(s) {
    filter(benchmark_dataset, signif == s, workflow == w) %>% .$Protein.IDs %>% as.character()
  })
})
sets <- list(x[[1]][[1]], x[[1]][[2]], x[[2]][[1]], x[[2]][[2]])
names(sets) <- c("pipeline_signif", "pipeline_no", "genedata_signif", "genedata_no")
areas <- lapply(sets, length) %>% unlist

venn_partitions <- get.venn.partitions(sets) %>% filter(..count.. > 0)  %>%
  select(pipeline_signif:genedata_no, ..count..) %>% rename(count = ..count..)

write.table(x = venn_partitions, file = file.path(output_dir, "venn_partitions.tsv"), col.names = T,row.names = F,quote = F,sep = "\t")

nono <- venn_partitions[1,]$..count..
yesno <- venn_partitions[2,]$..count..
noyes <- venn_partitions[4,]$..count..
yesyes <- venn_partitions[5,]$..count..
png(filename = file.path(output_dir, "venn.png"))
grid.draw(VennDiagram::draw.quad.venn(area1 = areas[1], area2 = areas[2], area3 = areas[3], area4 = areas[4],
                            n13 = yesyes, n24 = nono, n14 = yesno, n23 = noyes, n1234 = 0,n12 = 0,n34 = 0,n123 = 0,n124 = 0,n134 = 0,n234 = 0,
                            fill =  c("darkred", "red", "darkblue", "blue"), category = names(sets)))

dev.off()
             


link_pipeline_genedata <- function(genedata, pipeline) {
  match_row <- lapply(genedata$Protein.IDs, function(x) grep(pattern = x, x = pipeline$Protein.IDs))
  match_row <- lapply(match_row, function(x) ifelse(length(x) != 1, NA, x)) %>% unlist
  return(match_row)
}
 
genedata$pipeline_id <- link_pipeline_genedata(genedata, pipeline)
# Set protein ids from pipeline to the value of the matched one in genedata
# This way, the protein id from the protein group in pipeline that had a match in genedata is used
# pipeline$Protein.IDs[na.omit(genedata$pipeline)] <- genedata %>% filter(!is.na(pipeline)) %>% .$Protein.IDs




genedata <- genedata %>% filter(!is.na(pipeline_id))

pipeline <- pipeline[genedata$pipeline_id,]

estimate_benchmark <- data.frame(Protein.IDs = pipeline$Protein.IDs, pipeline = pipeline$estimate, genedata = genedata$log2FC)
estimate_benchmark <- estimate_benchmark[!(is.na(estimate_benchmark$pipeline) | is.na(estimate_benchmark$genedata)),]

# %>%  gather(workflow, estimate, -Protein.IDs)

correlation_plot <- ggplot(data = estimate_benchmark, aes(x = pipeline, y = genedata)) + geom_point() 
# bin_plot <- ggplot(data = estimate_benchmark, aes(x = pipeline, y = genedata)) + geom_bin2d()
correlation_plot

  # estimate_benchmark <- estimate_benchmark %>% filter(abs(pipeline) > 0.1)

linear_model <- lm(formula = "pipeline ~ genedata", data = estimate_benchmark, model = T)

summary(linear_model)

corr_eqn <- function(x,y, method="pearson", digits = 2) {
  corr_coef <- round(cor(x, y,method = method), digits = digits)
  letters <- c("rho")
  names(letters) <- "pearson"
  paste0("italic(", letters[method], ") == ", corr_coef)
}

correlation <- corr_eqn(estimate_benchmark$pipeline, estimate_benchmark$genedata)
correlation

correlation_plot + geom_text(x=-1, y=2, label=correlation, parse = T)
ggsave(filename = file.path(output_dir, "correlation_plot.png"))
