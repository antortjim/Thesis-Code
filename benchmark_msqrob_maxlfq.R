library("dplyr")
library("ggplot2")
setwd("MEGA/Master/Thesis/MaxLFQ/")
RSqM_signif <- read.table("RSqM_signif_peptides.txt",sep="\t", header=T) %>% filter(!is.na(taxonomy))
lfq_signif <- read.table("lfq_signif.tsv",sep="\t", header=T)
colnames(lfq_signif)[colnames(lfq_signif) == "taxonomy"] <- "taxonomy"

group_by(RSqM_signif, taxonomy) %>% summarise(count = n())
group_by(lfq_signif, taxonomy) %>% summarise(count = n())


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

library(patchwork)

p4+q4

