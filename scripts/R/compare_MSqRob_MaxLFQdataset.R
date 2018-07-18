library("readxl")
library("dplyr")
library("ggplot2")
library("cowplot")
library(pROC)
library(latex2exp)
library("viridis")
library("stringr")
library("tidyr")
library(VennDiagram)
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj",
                   ifelse(Sys.info()["user"] == "aoj", "/z/home/aoj",
                          "/home/antortjim/"))


data_dir <- ifelse(Sys.info()["user"] == "aoj", "thesis/genedata/maxlfq/",
                   "MEGA/Master/Thesis/Code/scripts/data")


report_dir <- file.path(home_dir, "MEGA/Master/Thesis/Report/plots/")



### 1. LOAD MQ+LFQ PROCESSED DATA
# Process published supplementary material
maxlfq_results <- read_xlsx(file.path(home_dir, data_dir, "mcp.M113.031591-1.xlsx"))
maxlfq_results <- maxlfq_results[,c(1:6, 33, grep(pattern = "Majority protein IDs", x = colnames(maxlfq_results)))]
lfq_intensities <- log2(maxlfq_results[,1:6] %>% as.matrix)
lfq_signif_list <- apply(lfq_intensities, 1, function(row) {
  x <- row[1:3]
  y <- row[4:6]
  t.test(x[!is.infinite(x)], y[!is.infinite(y)], alternative = "two.sided")
  })

lfq_signif_list <- lapply(lfq_signif_list, function(x) c("p.value" = x$p.value,
                                                         "estimate" = (x$estimate["mean of x"] - x$estimate["mean of y"]),
                                                         "method" = x$method,
                                                         "df" = x$parameter["df"]))

lfq_signif <- do.call(rbind, lfq_signif_list)

lfq_signif <- lfq_signif %>% as.data.frame
lfq_signif$taxonomy <- maxlfq_results$Taxonomy
lfq_signif$prot <- maxlfq_results$`Majority protein IDs`

colnames(lfq_signif)[2] <- "estimate"
lfq_signif$estimate <- lfq_signif$estimate %>% as.character %>% as.numeric
lfq_signif$p.value <- lfq_signif$p.value %>% as.character %>% as.numeric

lfq_signif <- lfq_signif %>% arrange(p.value)
lfq_signif$index <- 1:nrow(lfq_signif)
# lfq_signif$p.value <= (0.05 * (lfq_signif$index / nrow(lfq_signif)))
lfq_signif$qval <- p.adjust(p = lfq_signif$p.value, method = "fdr")
colnames(lfq_signif)[colnames(lfq_signif) == "Taxonomy"] <- "taxonomy"
group_by(lfq_signif, taxonomy) %>% summarise(count = n())

# Export to a standard format
#   write.table(x = lfq_signif, file = "data/lfq_signif.tsv", quote = F,sep = "\t",row.names = F,col.names = T)

### 2. LOAD MQ+MSqRob PROCESSED DATA
# Read MSqRob output
RSqM_signif <- read.table(file.path(home_dir, data_dir, "../../model/data/RSqM_signif_peptides.txt"),
                          sep="\t", header=T) %>% filter(!is.na(taxonomy))

# MQ_histogram_plot <- ggplot(data = lfq_signif, aes(x=estimate, fill = taxonomy)) + geom_histogram(position="dodge")
# MQ_density_plot <- ggplot(data = lfq_signif, aes(x=estimate, fill = taxonomy)) + geom_density(alpha=0.5)
# 
# MSqRob_histogram_plot <- ggplot(data = RSqM_signif, aes(x=estimate, fill = taxonomy)) + geom_histogram(position="dodge")
# MSqRob_density_plot <- ggplot(data = RSqM_signif, aes(x=estimate, fill = taxonomy)) + geom_density(alpha=0.5)

# plots <- list(MQ_histogram_plot, MQ_density_plot, MSqRob_histogram_plot, MSqRob_density_plot)
# plot_legend <- plots[[1]] %>% get_legend()
# plots <- plots %>% lapply(function(x) x + guides(fill=F) + scale_x_continuous(limits = c(-3,6)))
# plot_grid(plot_legend, plot_grid(plotlist = plots, nrow = 2, labels = c("MaxQuant", "", "MSqRob","")), nrow=2, rel_heights = c(1,5))


# RSqM_signif %>% filter(taxonomy != "Homo sapiens") %>% arrange(estimate) %>%
#   filter(estimate < 1)
# 
# RSqM_signif <- RSqM_signif %>% filter(se < 0.5)

# ggplot(RSqM_signif %>% filter(taxonomy == "Homo sapiens"), aes(x=estimate, y=se)) +
#   geom_point() + coord_cartesian(ylim=c(0,2))

# 
# semi_volcano_msqrob <- ggplot(RSqM_signif, aes(y=-log10(qval), x = abs(estimate), col = taxonomy)) +
#   geom_point(size=0.2) + geom_hline(yintercept=-log10(0.05))
# jitter_msqrob <-ggplot(RSqM_signif, aes(y=-log10(qval), x = taxonomy, group = taxonomy)) +
#   geom_jitter()
# volcano_msqrob <- ggplot(RSqM_signif, aes(y=-log10(qval), x = estimate, col = taxonomy)) +
#   geom_point() + facet_wrap(~taxonomy) + geom_vline(xintercept = 0, linetype="dashed")
# histogram_msqrob <- ggplot(RSqM_signif, aes(x = estimate, fill = taxonomy)) +
#   geom_histogram(position="dodge") + guides(fill=F)
# 
# semi_volcano_lfq <- ggplot(lfq_signif, aes(y=-log10(p.value), x = abs(estimate), col = taxonomy)) +
#   geom_point(size=0.2) + geom_hline(yintercept=-log10(0.05))
# 
# jitter_lfq <- ggplot(lfq_signif, aes(y=-log10(p.value), x = taxonomy, group = taxonomy)) +
#   geom_jitter()
# 
# volcano_lfq <- ggplot(lfq_signif, aes(y=-log10(p.value), x = estimate, col = taxonomy)) +
#   geom_point() + facet_wrap(~taxonomy) + geom_vline(xintercept = 0, linetype="dashed")
# 
# histogram_lfq <- ggplot(lfq_signif, aes(x = estimate, fill = taxonomy)) +
#   geom_histogram(position="dodge") + guides(fill=F)

# theme_set(theme_cowplot(font_size=12)) # reduce default font size
# plot_list <-list(semi_volcano_msqrob, semi_volcano_lfq) %>% lapply(function(x) x + guides(col=F))
# plot_grid(plot_list[[1]], plot_list[[2]],
#          labels = "AUTO")
 
# plot(
#   cumsum(RSqM_signif %>% arrange(qval) %>% .$taxonomy == "Homo sapiens") / nrow(RSqM_signif %>% filter(taxonomy == "Homo sapiens")),
#   cumsum(RSqM_signif %>% arrange(qval) %>% .$taxonomy != "Homo sapiens") / nrow(RSqM_signif %>% filter(taxonomy != "Homo sapiens")),
#   type="l", col = "red")
# lines(
#   cumsum(lfq_signif %>% arrange(p.value) %>% .$taxonomy == "Homo sapiens") / nrow(lfq_signif %>% filter(taxonomy == "Homo sapiens")),
#   cumsum(lfq_signif %>% arrange(p.value) %>% .$taxonomy != "Homo sapiens") / nrow(lfq_signif %>% filter(taxonomy != "Homo sapiens")),
#   type="l", col = "green")
# abline(a = 0, b = 1)
# 
# 
# plot(
#   cumsum(RSqM_signif %>% arrange(-estimate) %>% .$taxonomy == "Homo sapiens") / nrow(RSqM_signif_taxon %>% filter(taxonomy == "Homo sapiens")),
#   cumsum(RSqM_signif %>% arrange(-estimate) %>% .$taxonomy != "Homo sapiens") / nrow(RSqM_signif_taxon %>% filter(taxonomy != "Homo sapiens")),
#   type="l", col = "red")
# lines(
#   cumsum(lfq_signif %>% arrange(-estimate) %>% .$taxonomy == "Homo sapiens") / nrow(lfq_signif %>% filter(taxonomy == "Homo sapiens")),
#   cumsum(lfq_signif %>% arrange(-estimate) %>% .$taxonomy != "Homo sapiens") / nrow(lfq_signif %>% filter(taxonomy != "Homo sapiens")),
#   type="l", col = "green")
# abline(a = 0, b = 1)


### 3. LOAD Compomics+MSqRob PROCESSED DATA
RSqM_signif_compomics <- read.table(file.path(home_dir, data_dir, paste0("RSqM_signif_taxon", ".tsv")),
                    sep = "\t", header = T, stringsAsFactors = F)

# source(file.path(home_dir, data_dir, "../R/check_organism.R"))
# organism_check <- check_organism(RSqM_signif_compomics$Protein.IDs, split = ", ")
# 
# RSqM_signif_compomics <- left_join(RSqM_signif_compomics,
#          data.frame(Protein.IDs = RSqM_signif_compomics$Protein.IDs, taxonomy = organism_check[[2]]),
#          by = "Protein.IDs") %>%
#  filter(!is.na(taxonomy))

RSqM_signif_taxon <- RSqM_signif_compomics %>% filter(abs(estimate) > 1e-6)
RSqM_signif_taxon <- RSqM_signif_compomics 

# write.table(x = RSqM_signif_taxon, file = file.path(export_folder, paste0("RSqM_signif_taxon", suffix, ".tsv")),
#             sep = "\t", col.names = T, row.names = F, quote = F)
 

### 4. PLOT DATA
nbins <- 100
# xaxis <- TeX('$|\\log_2(FC)|$')
# yaxis <- TeX('$-\\log_{10}(qval)$')
# 
# theme_set(theme_gray())
# volcano_compomics <- ggplot(
#   data = RSqM_signif_taxon,
#   mapping = aes(y=-log10(qval), x = estimate, col = taxonomy)) +
#   geom_point() +
#   geom_hline(yintercept=-log10(0.05)) +
#   xlab(xaxis) + ylab(yaxis) + guides(col=guide_legend(title="Species"))
# 
# semi_volcano_compomics <- ggplot(
#   data = RSqM_signif_taxon,
#   mapping = aes(y=-log10(qval), x = abs(estimate), col = taxonomy)) +
#   geom_point() +
#   geom_hline(yintercept=-log10(0.05)) +
#   scale_x_continuous(breaks=0:8) +
#   xlab(xaxis) + ylab(yaxis) + guides(col=guide_legend(title="Species"))
# 
# 
# histogram_compomics <- ggplot(
#   data = RSqM_signif_taxon,
#   mapping = aes(x =estimate, fill = taxonomy)) +
#   geom_histogram(bins=nbins, position="dodge") +
#   guides(fill=guide_legend(title="Species"))
# 
# volcano_compomics
# semi_volcano_compomics
# histogram_compomics


# ggplot(data = RSqM_signif_taxon, aes(x=estimate, fill=taxonomy)) +
#   geom_histogram(bins=50, position="dodge") +
#   coord_cartesian(xlim = c(-5,5)) +
#   scale_x_continuous(breaks = -5:5)

# ggsave(file.path(report_dir, "histogram_pipeline.png") , height = 7, width = 7)


pipelines <- c("Comp+Rob", "MQ+LFQ", "MQ+Rob"
               # , "MQ+MSBay"
               )

# MSBayQ <- read.table(file = file.path(data_dir, "../../model/data/MSBayQ.tsv"), header=T, sep = "\t")
# colnames(MSBayQ)[1] <- "protein"
# ROPE <- 0.4
# MSBayQ$significance <- ifelse((MSBayQ$hpd_2.5 > ROPE), 0, 1)


tools_data <- rbind(
  cbind(select(RSqM_signif_taxon, estimate, qval, taxonomy, Protein.IDs) %>% rename(prot = Protein.IDs, significance = qval), method = pipelines[1]),
  cbind(select(lfq_signif, estimate, qval, taxonomy, prot) %>% rename(significance=qval), method = pipelines[2]),
  # cbind(select(MSBayQ, mean, significance, Organism) %>% rename(estimate=mean, taxonomy=Organism), method = pipelines[4]),
  cbind(select(RSqM_signif, estimate, qval, taxonomy, protein) %>% rename(prot = protein, significance = qval), method = pipelines[3])
  )

tools_data$taxonomy <- tools_data$taxonomy %>% as.character()
tools_data[tools_data$taxonomy == "Escherichia coli (strain K12)", "taxonomy"] <- "E. coli"

colnames(tools_data) <- c("log2FC", "signif", "Organism", "prot", "Tool")


make_plot_data <- function(curve_data, curve="roc") {
  control_data <- data.frame(x = NULL, y = NULL, method = NULL)
  pipelines <- curve_data$Tool %>% unique
  aucs <- rep(0, length(pipelines))
  i <- 1
  for(pip in pipelines) {
    print(pip)
    TN <- (curve_data %>% filter(Tool == pip) %>% arrange(signif) %>% .$Organism) == "Homo sapiens"
    TP <- (curve_data %>% filter(Tool == pip) %>% arrange(signif) %>% .$Organism) != "Homo sapiens"
    
    if(curve=="roc") {
      control_data <- rbind(control_data,
                             data.frame(x=(cumsum(TN) / sum(TN)), y=(cumsum(TP) / sum(TP)), Tool = pip)
      )
      control_data <- rbind(control_data, data.frame(x=1,y=1, Tool=pip))
      
      aucs[i] <- roc(curve_data %>% filter(Tool == pip) %>% .$Organism == "Homo sapiens",
                     curve_data %>% filter(Tool == pip) %>% .$signif
      )$auc %>% round(2)
      
    } else if (curve=="prc") {
      control_data <- rbind(control_data,
                             data.frame(
                               x=cumsum(TP)/sum(TP),
                               y = cumsum(TP)/(cumsum(TP) + cumsum(TN)),
                               Tool = pip
                               )
                            )
      control_data <- rbind(control_data, data.frame(x=1,y=0, Tool=pip))
    }
    
    i <- i + 1
  }
  return(list(control_data = control_data, aucs = aucs))
}

plot_curve <- function(curve_data, curve="roc") {
  
  fpr <- TeX("False Positive Rate ($\\frac{FP}{FP + TN}$)")
  fnr <- TeX("True Positive Rate ($\\frac{TP}{FN + TP}$)")
  precision <- TeX("Precision ($\\frac{TP}{TP + FP}$)")
  
  xcoord <- Inf
  ycoord <- -Inf
  pipelines <- curve_data$Tool %>% unique
  n <- 4
  
  palette <- c("#b42f32", "#e3e3cd", "#eba631", "black")[1:length(pipelines)]
  
  dd <- make_plot_data(curve_data = curve_data, curve=curve) 
  curve_data <- dd[[1]]
  aucs <- dd[[2]]
  
  curve_plot <- ggplot(data = curve_data, aes(x=x,y=y, col=Tool)) +
    geom_line(size=1.2) +
    theme_bw() +
    theme(legend.position = "top") +
    scale_color_manual(values = palette) +
    coord_cartesian(ylim=c(0,1))
  
  if(curve == "roc") {
    for (i in 1:length(pipelines)) {
      annot <- str_pad(string = aucs[i], width = 4, side = "right", pad = "0")
      
      curve_plot <- curve_plot +
      annotate(geom = "label", x = xcoord, y = ycoord, label=paste0("AUC: ", annot),
               fill=palette[i], color="white", fontface="bold",
               hjust = 1.2, vjust = -1*i/2 -(i-0.5))
    }
    curve_plot <- curve_plot + xlab("FPR") + ylab("TPR") + geom_abline(slope=1, intercept=0, linetype="dashed")
  } else {
    curve_plot <- curve_plot + xlab("TPR") + ylab("Precision") + geom_abline(slope=-1, intercept=1, linetype="dashed")
  }
  curve_plot <- curve_plot + theme_bw(base_size=15)
  return(curve_plot)
}
ROC_curve <- plot_curve(tools_data)
ROC_curve_FC <- plot_curve(tools_data %>% filter(log2FC > 1))
PRC_curve <- plot_curve(tools_data, curve = "prc")
PRC_curve_FC <- plot_curve(tools_data %>% filter(log2FC > 1), curve = "prc")
curves <- list(ROC_curve, ROC_curve_FC, PRC_curve, PRC_curve_FC) %>%
  lapply(., function(x) x + guides(col=F))

curves_plot <- plot_grid(
  get_legend(ROC_curve + theme(legend.direction = "horizontal", legend.title = element_text(size=20), legend.text = element_text(size=15))),
          plot_grid(plotlist = curves, nrow=2,ncol=2, labels="AUTO"), nrow=2, rel_heights = c(1,10))

curves_plot
ggsave(filename = file.path(report_dir, "curves_plot.png"), plot = curves_plot, height=8, width=8)


tools_data <- tools_data %>% filter(!is.na(signif))
tools_data %>% group_by(Tool, Organism, signif < 0.05)  %>% summarise(count = n())
tools_data$Organism <- factor(tools_data$Organism, levels=c("Homo sapiens", "E. coli"))

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
palette <- rev(gg_color_hue(2))

ylab_tex <- TeX("$-log_{10}$(\\textit{q})")
plot(ylab_tex)
box_plot <- ggplot(tools_data, aes(fill = Organism, x = Organism, y = -log10(signif))) +
  # geom_jitter() +
  geom_boxplot() +
  guides(fill=F) +
  facet_wrap(~Tool) +
  labs(y = ylab_tex) +
  scale_fill_manual(values = palette)

theme_set(theme_bw(base_size = 15) + theme(legend.title = element_text(size = 20)))
histogram_plot <- ggplot(tools_data, aes(fill = Organism, x = log2FC)) +
  # geom_jitter() +
  geom_histogram(position="dodge", bins=nbins) +
  facet_wrap(~Tool) +
  theme(legend.position = "top") +
  coord_cartesian(xlim=c(-3,3), ylim=c(0,750)) +
  scale_x_continuous(breaks=(-3:3)) +
  labs(y = "Frequency")+
  scale_fill_manual(values = palette)

compomics_histogram_plot <- ggplot(filter(tools_data, Tool == "Compomics+MSqRob"), aes(fill = Organism, x = log2FC)) +
  # geom_jitter() +
  geom_histogram(position="dodge", bins=nbins) +
  # facet_wrap(~Tool) +
  theme(legend.position = "top") +
  coord_cartesian(xlim=c(-3,3), ylim=c(0,750)) +
  scale_x_continuous(breaks=(-3:3)) +
  labs(y = "Frequency")+
  scale_fill_manual(values = palette) +
  theme(legend.position = "top")


# ggplot(MSBayQ %>% filter(n_peptides>3), aes(fill = Organism, x = mean)) +
#   # geom_jitter() +
#   geom_histogram(position="dodge", bins=60) +
#   theme(legend.position = "top") +
#   # coord_cartesian(xlim=c(-3,3)) +
#   # scale_x_continuous(breaks=(-3:3)) +
#   labs(y = "Frequency")+
#   scale_fill_manual(values = palette)



volcano_plot <- ggplot(tools_data, aes(col = Organism, x = log2FC, y = -log10(signif))) +
  # geom_jitter() +
  geom_point(size=0.2) +
  facet_wrap(~Tool) +
  theme(legend.position = "top") +
  coord_cartesian(xlim=c(-3,3)) +
  scale_x_continuous(breaks=(-3:3)) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  labs(y = ylab_tex) +
  scale_color_manual(values = palette)


combined_plot <- plot_grid(histogram_plot,
          volcano_plot + theme(legend.position = "none") + geom_hline(yintercept = -log10(0.05), linetype="dashed"),
          box_plot + geom_hline(yintercept = -log10(0.05), linetype="dashed"),
          nrow=3, labels="AUTO", align = "v")

combined_plot

ggsave(filename = file.path(report_dir, "combined_plot.png"), plot = combined_plot, height=10, width=8)
ggsave(filename = file.path(report_dir, "histogram_compomics.png"), plot = compomics_histogram_plot, height=5, width=10)


true_data <- tools_data %>% filter(log2FC > 1, signif < 0.05, Organism == "E. coli") %>% select(prot, Tool)
strsplit(true_data$prot[grep(";", true_data$prot)], split = ";") %>%
  lapply(function(x) strsplit(x, split = "-") %>% lapply(function(y) {y[1]}) %>%
           unlist %>% unique %>% paste(., collapse = ";")) %>% unlist

true_positives <- lapply(tools_data$Tool %>% unique, function(tool) true_data %>% filter(Tool == tool) %>% .$prot)

shared_positives1 <- true_positives[[1]][true_positives[[1]] %in% true_positives[[2]]]
shared_positives2 <- true_positives[[1]][true_positives[[1]] %in% true_positives[[3]]]
shared_positives <- shared_positives1[shared_positives1 %in% shared_positives2]

for (tool in pipelines) {
  write.table(x = ,
              file = file.path(home_dir, data_dir, paste0(tool, "_true.txt")),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

pipelines

v <- venn.diagram(x = lapply(pipelines, function(x) true_data[true_data$Tool ==x, "prot"]),
                  filename = file.path(report_dir, "vennDiagram.png"), imagetype = "png",
                  # filename = NULL,
                  category.names = tools_data$Tool %>% unique,
                          fill = c("#ff0000",
                                   # "#5b00ff",
                                   "#f5ff6e",
                                   "#00ff0c"),
                  alpha = 0.5, lwd=0, cex = 1.7, cat.cex=1.3,
                  cat.just=list(c(0,6) , c(1.4,6) , c(.5,-6)))

grid.newpage() & grid.draw(v)

filter(tools_data, prot == "P02919")
shared_data <- tools_data %>% filter(prot %in% shared_positives) %>% select(log2FC, prot, Tool)

pipeline_pairs <- combn(pipelines, 2)
corr_plot_data <- data.frame(prot=NULL, x=NULL, y=NULL, pair=NULL)
for (i in 1:ncol(pipeline_pairs)) {
  current_pair <- pipeline_pairs[,i]
  xx <- spread(filter(shared_data, Tool %in% current_pair), Tool, log2FC)
  colnames(xx)[2:3] <- c("x", "y")
  xx$Tool <- paste(current_pair, collapse = " - ")
  corr_plot_data <- rbind(corr_plot_data, xx)
}


correlation <- group_by(corr_plot_data, Tool) %>% summarise(corr = cor(x, y))

xcoord <- Inf
ycoord <- -Inf

rho <- sapply(round(correlation$corr, digits=2), function(co) {
  paste0("r = ", co)
  }) %>% as.list %>% unlist


pipeline_pairs_char <- apply(pipeline_pairs, 2, function(x) paste(x, collapse = " - "))

annot_data <- data.frame(Tool = pipeline_pairs_char,
                         rho = rho)

facet_labels <- pipeline_pairs_char
names(facet_labels) <- 1:3

correlation_plot <- ggplot(data = corr_plot_data, aes(x=x,
                                                      # col=Tool,
                                                      y=y)) +
  geom_point() +
  facet_wrap(~Tool
             # labeller = labeller(Tool = facet_labels)
             ) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  coord_cartesian(xlim = 0.5:4, ylim = 0.5:4) +
  labs(x = "First pipeline", y = "Second pipeline")


correlation_plot <- correlation_plot +
           geom_text(data = annot_data,
                     aes(x = xcoord, y = ycoord, col=NULL, label=rho),
                     # label="hola",
  # fill=palette[i], color="white", fontface="bold",
  hjust = 1.6, vjust =-2.2, size=7) +
  theme_bw(base_size=20) +
  theme(legend.position = "top")
correlation_plot

ggsave(filename = file.path(report_dir, "correlation_plot.png"),
       plot = correlation_plot, height=6, width=8)


RSqM_signif_taxon %>% group_by(Organism) %>% summarise(sum((log2FC > log2(2.5) & qval <  0.05), na.rm=T))


RSqM_table <- RSqM_signif_taxon[complete.cases(RSqM_signif_taxon),]
RSqM_table <- RSqM_table[(grep(x = RSqM_table$Protein.IDs, pattern = ", ", invert = T)),] %>%
  arrange(qval) %>%
  select(estimate, qval, Protein.IDs, taxonomy) %>%
  rename(log2FC = estimate, Protein=Protein.IDs, Organism=taxonomy)
rownames(RSqM_table) <- NULL
RSqM_table$Organism <- as.character(RSqM_table$Organism)
RSqM_table$Organism[RSqM_table$Organism == "Escherichia coli (strain K12)"] <- "E. coli"
RSqM_table <- RSqM_table %>% head(10)

# https://stackoverflow.com/questions/16579562/r-knitr-xtable-alternating-row-colors
rws <- seq(1, (nrow(RSqM_table)-1), by = 2)
col <- rep("\\rowcolor[gray]{0.95}", length(rws))
print(xtable::xtable(RSqM_table, digits=-2,
                     align = c(rep("l", 4), ">{\\itshape}l")), booktabs = TRUE,
               add.to.row = list(pos = as.list(rws), command = col))

      # add.to.col = list(pos = list(colnames(RSqM_table) == "Organism"), command = "\\textit"))


