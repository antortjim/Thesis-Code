library("readxl")
library("dplyr")
library("ggplot2")
library("cowplot")
setwd("~/MEGA/Master/Thesis/Code/")

# Process published supplementary material
maxlfq_results <- read_xlsx("MaxLFQ/mcp.M113.031591-1.xlsx")
maxlfq_results <- maxlfq_results[,c(1:6, 33)]
lfq_intensities <- maxlfq_results[,1:6] %>% as.matrix
lfq_signif_list <- apply(lfq_intensities, 1, function(x) t.test(x[1:3], x[4:6], alternative = "two.sided"))
lfq_signif_list <- lapply(lfq_signif_list, function(x) c("p.value" = x$p.value,
                                                         "estimate" = log2(x$estimate["mean of x"] / x$estimate["mean of y"]),
                                                         "method" = x$method,
                                                         "df" = x$parameter["df"]))
lfq_signif <- do.call(rbind, lfq_signif_list)

lfq_signif <- lfq_signif %>% as.data.frame
lfq_signif$taxonomy <- maxlfq_results$Taxonomy

colnames(lfq_signif)[2] <- "estimate"
lfq_signif$estimate <- lfq_signif$estimate %>% as.character %>% as.numeric
lfq_signif$p.value <- lfq_signif$p.value %>% as.character %>% as.numeric

lfq_signif <- lfq_signif %>% arrange(p.value)
lfq_signif$index <- 1:nrow(lfq_signif)
lfq_signif$p.value <= (0.05 * (lfq_signif$index / nrow(lfq_signif)))
lfq_signif$qval <- p.adjust(p = lfq_signif$p.value, method = "fdr", n = length(lfq_signif$p.value))

colnames(lfq_signif)[colnames(lfq_signif) == "Taxonomy"] <- "taxonomy"
write.table(x = lfq_signif, file = "data/lfq_signif.tsv", quote = F,sep = "\t",row.names = F,col.names = T)
# Export to a standard format

# Read MSqRob output
RSqM_signif <- read.table("data/RSqM_signif_peptides.txt",sep="\t", header=T) %>% filter(!is.na(taxonomy))

group_by(RSqM_signif, taxonomy) %>% summarise(count = n())
group_by(lfq_signif, taxonomy) %>% summarise(count = n())

# Filter to keep proteins exactly the same

MQ_histogram_plot <- ggplot(data = lfq_signif, aes(x=estimate, fill = taxonomy)) + geom_histogram(position="dodge")
MQ_density_plot <- ggplot(data = lfq_signif, aes(x=estimate, fill = taxonomy)) + geom_density(alpha=0.5)

MSqRob_histogram_plot <- ggplot(data = RSqM_signif, aes(x=estimate, fill = taxonomy)) + geom_histogram(position="dodge")
MSqRob_density_plot <- ggplot(data = RSqM_signif, aes(x=estimate, fill = taxonomy)) + geom_density(alpha=0.5)

plots <- list(MQ_histogram_plot, MQ_density_plot, MSqRob_histogram_plot, MSqRob_density_plot)
plot_legend <- plots[[1]] %>% get_legend()
plots <- plots %>% lapply(function(x) x + guides(fill=F) + scale_x_continuous(limits = c(-3,6)))
plot_grid(plot_legend, plot_grid(plotlist = plots, nrow = 2, labels = c("MQ", "", "MSqRob","")), nrow=2, rel_heights = c(1,5))


RSqM_signif %>% filter(taxonomy != "Homo sapiens") %>% arrange(estimate) %>%
  filter(estimate < 1)

p1 <- ggplot(RSqM_signif, aes(y=-log10(qval), x = abs(estimate), col = taxonomy)) +
  geom_point(size=0.2) + geom_hline(yintercept=-log10(0.05))
p2 <-ggplot(RSqM_signif, aes(y=-log10(qval), x = taxonomy, group = taxonomy)) +
  geom_jitter()
p3 <- ggplot(RSqM_signif, aes(y=-log10(qval), x = estimate, col = taxonomy)) +
  geom_point() + facet_wrap(~taxonomy) + geom_vline(xintercept = 0, linetype="dashed")
p4 <- ggplot(RSqM_signif, aes(x = estimate, fill = taxonomy)) +
  geom_histogram(position="dodge") + guides(fill=F)

q1 <- ggplot(lfq_signif, aes(y=-log10(qval), x = abs(estimate), col = taxonomy)) +
  geom_point(size=0.2) + geom_hline(yintercept=-log10(0.05))

q2 <- ggplot(lfq_signif, aes(y=-log10(qval), x = taxonomy, group = taxonomy)) +
  geom_jitter()

q3 <- ggplot(lfq_signif, aes(y=-log10(qval), x = estimate, col = taxonomy)) +
  geom_point() + facet_wrap(~taxonomy) + geom_vline(xintercept = 0, linetype="dashed")
q4 <- ggplot(lfq_signif, aes(x = estimate, fill = taxonomy)) +
  geom_histogram(position="dodge") + guides(fill=F)

p1
q1
