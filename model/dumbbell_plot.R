library(ggplot2); library(dplyr); library(tidyr); library(ggalt); library(cowplot)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
palette <- rev(gg_color_hue(2))


make_dumbbell_plot <- function(MSBayQ, n_pep=NULL, n=5) {
  
  print(n)
  # if (!is.null(n_pep)) MSBayQ <- MSBayQ[MSBayQ$n_peptides == n_pep,]
  MSBayQ <- MSBayQ %>%
    group_by(n_peptides, Organism) %>%
    filter(row_number() < 6)
  
  number_facets <- length(unique(MSBayQ$n_peptides))
  
  MSBayQ <- MSBayQ %>% arrange(n_peptides, Organism)
  MSBayQ$protein <- rep(seq(0,(n/2-.5), by=.5), times=number_facets)
      
  # selected <- lapply(selected, function(x) sample(x, min(n, length(x)))) %>% unlist
  
  # MSBayQ <- MSBayQ[MSBayQ$protein %in% selected, ]
  MSBayQ$protein <- as.character(MSBayQ$protein) %>% lapply(function(x) strsplit(x, split=";") %>% unlist %>% .[1] %>% substr(1,6)) %>% unlist
  #     print(MSBayQ))
  p <- ggplot(data=MSBayQ, aes(x=hpd_2.5, xend=hpd_97.5, y = protein, color=Organism)) +
    geom_vline(xintercept=-0.4, col=palette[1], linetype="dashed") +
    geom_vline(xintercept=0.4, col=palette[1], linetype="dashed") +
    geom_vline(xintercept=log2(3), col=palette[2], linetype="dashed") +
    geom_dumbbell(size=0.75,
                  point.colour.l="#0e668b") +
    geom_point(aes(x=mean, y = protein)) +
    labs(x="log2FC 95% HPDI", 
         y="Proteins")
  
  if(!is.null(n_pep)) {
  p <- p + ggtitle(paste0("Performance. ", n_pep, " peptides"))
}
 p <- p + theme_bw(base_size=15) + theme(plot.title = element_text(hjust=0.5, face="bold"),
                axis.text.y=element_blank(),
          # plot.background=element_rect(fill="#f7f7f7"),
          # panel.background=element_rect(fill="#f7f7f7"),
          # panel.grid.minor=element_blank(),
          # panel.grid.major.y=element_blank(),
          # panel.grid.major.x=element_line(),
          axis.ticks.y=element_blank(),
          legend.position="top") +
    scale_x_continuous(breaks=seq(from=-1, to=3, by=1), limits=c(-1,3)) +
    facet_wrap(~n_peptides, nrow = number_facets/2)
  return(p)
}

make_error_hdi_interval_plot <- function(MSBayQ, n_pep=NULL) {
  if (!is.null(n_pep)) MSBayQ <- MSBayQ[MSBayQ$n_peptides == n_pep,]
  MSBayQ$true <- 0
  MSBayQ[MSBayQ$Organism == "Escherichia coli (strain K12)","true"] <- log2(3)
  print(cor(x=MSBayQ$`hpd_97.5` - MSBayQ$`hpd_2.5`, y=abs(MSBayQ$mean-MSBayQ$true)))
  ggplot(data = MSBayQ, mapping = aes(x = (`hpd_97.5` - `hpd_2.5`), y = abs(mean - true), col = Organism)) + geom_point()
  
}


MSBayQ <- read.table("data/MSBayQ.tsv", header=T, sep = "\t")
colnames(MSBayQ)[1] <- "protein"
print(colnames(MSBayQ))

p1 <- make_dumbbell_plot(MSBayQ[MSBayQ$n_peptides %in% c(2,3,4,6,7,10),], n_pep = NULL, n=10) 
p1
ggsave(paste0("../../Report/plots/performance.eps"), height = 7, width=5, plot = p1)
ggsave(paste0("../../Report/plots/performance.png"), height = 7, width=5, plot = p1)
