library(ggplot2); library(dplyr); library(tidyr); library(ggalt); library(cowplot)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
palette <- rev(gg_color_hue(2))


make_dumbbell_plot <- function(BayesQuant, n_pep=NULL, n=5,sep=0.25) {
  
  # if (!is.null(n_pep)) BayesQuant <- BayesQuant[BayesQuant$n_peptides == n_pep,]
  bayesquant <- bayesquant %>%
    group_by(n_peptides, Organism) %>%
    filter(row_number() < n+1)
  
  number_facets <- length(unique(bayesquant$n_peptides))
  print(number_facets)
  print(n)
  
  bayesquant <- bayesquant %>% arrange(n_peptides, Organism)
  bayesquant$protein <- rep(seq(from = 0, by=sep, length.out = n*2), times=number_facets)
      
  # selected <- lapply(selected, function(x) sample(x, min(n, length(x)))) %>% unlist
  
  # BayesQuant <- BayesQuant[BayesQuant$protein %in% selected, ]
  bayesquant$protein <- as.character(bayesquant$protein) %>% lapply(function(x) strsplit(x, split=";") %>% unlist %>% .[1] %>% substr(1,6)) %>% unlist
  #     print(BayesQuant))
  p <- ggplot(data=bayesquant, aes(x=hpd_2.5, xend=hpd_97.5, y = protein, color=Organism)) +
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

make_error_hdi_interval_plot <- function(BayesQuant, n_pep=NULL) {
  if (!is.null(n_pep)) BayesQuant <- BayesQuant[BayesQuant$n_peptides == n_pep,]
  BayesQuant$true <- 0
  BayesQuant[BayesQuant$Organism == "Escherichia coli (strain K12)","true"] <- log2(3)
  print(cor(x=BayesQuant$`hpd_97.5` - BayesQuant$`hpd_2.5`, y=abs(BayesQuant$mean-BayesQuant$true)))
  ggplot(data = BayesQuant, mapping = aes(x = (`hpd_97.5` - `hpd_2.5`), y = abs(mean - true), col = Organism)) + geom_point()
  
}


BayesQuant <- read.table("data/bayesquant_res.tsv", header=T, sep = "\t")
colnames(BayesQuant)[1] <- "protein"
print(colnames(BayesQuant))

bayesquant <- BayesQuant[BayesQuant$n_peptides %in% c(2,10),]

p1 <- make_dumbbell_plot(bayesquant, n_pep = NULL, n=5,sep=0.05) 
p1
ggsave(paste0("../../Report/plots/performance.eps"), height = 7, width=5, plot = p1)
ggsave(paste0("../../Report/plots/performance.png"), height = 7, width=5, plot = p1)
