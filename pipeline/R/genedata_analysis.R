library("dplyr")
library("optparse")
library("ggplot2")
library("tidyr")

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="thp1"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name
output_dir <- opt$output_dir
output_dir <- ifelse(output_dir %>% is.null, file.path(root_dir, exp_name, "genedata_processed"), output_dir)

#######################################################################################
## Custom function
#######################################################################################

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


# pipeline_settings <- read.table(file = file.path(root_dir, exp_name, "pipeline_settings.txt"), sep = ":")

#######################################################################################
## Analysis parameters
#######################################################################################

dep_criteria <- c("BOTH")
palette <- rep("red", 4)


log.fold.change.threshold <- 1
p.value.threshold <- -log10(0.05)
correlation.threshold <- 0.96
#excluded_proteins <- "P33BJB"
analysis_conditions <- c("Experiment", "Fraction")
## Compute possible comparisons
comparisons <-matrix(c("pLPS_pellet", "mLPS_pellet"))
indices <- lapply(comparisons[,1], function(x) {
  grep(x, colnames(data_spread))
})


#######################################################################################
## Read data and save available conditions
#######################################################################################

data <- read.table(file = file.path(output_dir, "dataset.tsv"), stringsAsFactors = F, header = T, row.names = NULL)
data %>% group_by_(analysis_conditions, "Replicate") %>% summarise(count = n())

# #######################################################################################
# ##  COMPLETE
# #######################################################################################
# 
# counts_condition <- data %>%
#   group_by_(analysis_conditions, "Name") %>% summarise(count = sum(!is.na(Quantification)))
# 
# counts_condition
# 
# p <- ggplot(counts_condition ,
#        aes(x = Name, fill=Experiment, y = count)) + geom_bar(stat = "identity")
#   # scale_x_continuous(breaks = 0:max(counts_condition$count)) +
# p
# 
# ggsave2(filename = "Quantification_per_sample", plot = p)
# 
# # data <- data %>% filter(name %in% (data %>% group_by(condition, name) %>%
# #                                      summarise(count = n()) %>%
# #                                      filter(count > 4) %>% .$name))
# 
# # data[is.na(data$Quantification),"Quantification"] <- 0
# 
# #######################################################################################
# ## Boxplots of the protein Quantification distributions across samples
# #######################################################################################
# 
# p <- ggplot(data,
#        aes(x = Replicate, y = Quantification, fill=Experiment, group=Replicate)) +
#   geom_boxplot() + facet_wrap(~Experiment, nrow = 1)
# p
# 
# ## ADD NUMBER OF PROTEINS FOR EACH SAMPLE
# 
# ggsave2(filename = "Quantification_boxplots.png", plot = p)
# 
# 
# 
# # p <- ggplot(data %>% filter(!(Name %in% excluded_proteins)), aes(x = Experiment, y = log(Quantification))) +
# #   geom_boxplot()
# # p
# # ggsave2(filename = "Quantifications_condition_boxplot.png", plot = p)
# 
# p <- ggplot(data %>% filter(!(Protein.ID %in% excluded_proteins)), aes(x = Quantification, col = Experiment)) +
#   geom_density(alpha = 0.5) + guides(col=F)
# p
# ggsave2(filename = "Quantification_density.png", plot = p)
# 
# 
# 
# #######################################################################################
# ## Scatterplot of the proteins Quantification across conditions (one datapoint per replicate)
# #######################################################################################
# 
# ggplot(data %>% filter(!(name %in% excluded_proteins)), aes(x = name, y = Quantification, col=Experiment)) +
#   geom_point() +
#   facet_wrap(~condition) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# #######################################################################################
# ## Scatter plot of mean and coefficient of variation
# #######################################################################################
# 
# cv_data <- data %>%
#   group_by(condition, name) %>%
#   summarise(cv = sd(Quantification, na.rm = T) / mean(Quantification, na.rm = T),
#             mu = mean(Quantification, na.rm=T),
#             sigma = sd(Quantification, na.rm=T))
# 
# ggplot(cv_data, aes(x = mu, y = cv)) + geom_point() + facet_wrap(~condition)
# 
# # ggplot(cv_data %>% arrange(condition, -mu), aes(x = name, y = mu)) + geom_bar(stat="identity") + geom_errorbar(aes(x = name,
# #                                                                     ymin = mu - sigma,
# #                                                                     ymax = mu + sigma)) +
# #   facet_wrap(~condition)
# # 
# # data$name %>% as.character %>% unique
# 
# #######################################################################################
# ## Proteins quantified per sample
# #######################################################################################
# 
# quant_count <- data %>% na.omit %>% group_by(file) %>%
#   summarise(count = n()) %>% mutate(file = as.factor(file))
# ggplot(data = quant_count, aes(x = file, y = count)) + geom_bar(stat = "identity") + labs(y = "Protein count")
# 
# 
# # data %>% filter(file == 1) %>% .$Quantification
# # 
# # join_data <- inner_join(data[data$file == 1, c("name", "Quantification")], data[data$file == 2, c("name", "Quantification")], by = "name")
# # join_data[join_data %>% is.na] <- 0.1
# # 
# # M <- log2(join_data$Quantification.x) - log2(join_data$Quantification.y)
# # A <- 0.5*(log2(join_data$Quantification.x) + log2(join_data$Quantification.y))
# # 
# # plot(A, M, pch = 16)


# #######################################################################################
# ## Protein quantities across conditions
# #######################################################################################
# 
average_dat <- data %>%
  group_by(Experiment, Protein.ID) %>% summarise(average = mean(Quantification)) %>% arrange(Protein.ID)

# Keep useful features 
data_spread <- data %>% select(Experiment, Replicate, Quantification, Protein.ID)
# Keep only pellet
data_spread <- data_spread[data_spread$Experiment %>% grep(pattern = "pellet"),]
# Create a new sample field to identify each individual sample
data_spread$Sample <- paste0(data_spread$Experiment, "_", data_spread$Replicate)
data_spread <- data_spread %>% select(-Experiment, -Replicate)
# Make data spread
data_spread <- data_spread %>% spread(Sample, Quantification)
# Make protein.id rownames instead of separate field
rownames(data_spread) <- data_spread$Protein.ID
data_spread <- data_spread %>% select(-Protein.ID)
# Remove proteins seen in only one sample or only one condition

data_spread <- data_spread[!(data_spread %>% is.na %>% rowSums() >= 5) | (data_spread[,indices[[1]]] %>% is.na %>% rowSums == 0) | (data_spread[,indices[[2]]] %>% is.na %>% rowSums) ,]

data_spread_bad <- data_spread[!((data_spread[,indices[[1]]] %>% is.na %>% rowSums) <= 1 & (data_spread[,indices[[2]]] %>% is.na %>% rowSums) <= 1),]
## Implement method to work with proteins seen in both conditions but only once in at least one of them

data_spread <- data_spread[((data_spread[,indices[[1]]] %>% is.na %>% rowSums) <= 1 & (data_spread[,indices[[2]]] %>% is.na %>% rowSums) <= 1),]


impute_quants <- function(data_spread, indices=list(1:3, 4:6)) {
  
  for (j in 1:length(indices)) {
    missing_mean <- rowMeans(data_spread[,indices[[j]]], na.rm = T)
    missing_data <- which(data_spread %>% is.na, arr.ind = T)
    missing_data <- missing_data[missing_data[,2] >= indices[[j]][1] & missing_data[,2] <= tail(indices[[j]],1),]
    imputed <- apply(missing_data, 1, function(x) max(0, rnorm(n = 1, mean = missing_mean[x[1]], sd = sd(data_spread[x[1],indices[[j]]], na.rm = T))))
    for(i in 1:nrow(missing_data)) {
      data_spread[missing_data[i,1], missing_data[i,2]] <- imputed[i]
    }
  }
  
  return(data_spread)
}

data_spread <- impute_quants(data_spread)


## Keep proteins that are quantified in all samples
protein_names <- rownames(data_spread)



## Compute fold change and significance
# fc <- 1:ncol(comparisons) %>% lapply(function(i) {
#   mean(data_spread[, comparisons[1, i]]) / mean(data_spread[, comparisons[2, i]])
#   })
p_values <- matrix(ncol = ncol(comparisons), nrow = length(protein_names))
t_tests <- list()

for (j in 1:ncol(comparisons)) {
  i <- 1
  for(protein in protein_names) {
    #print(protein)
    t_test <- tryCatch({
      t.test(data_spread[i,indices[[1]]], data_spread[i,indices[[2]]], alternative = "two.sided")
      },
      error = function(e) {
        NA
    })
    # p_values[i, j] <- p_value
    t_tests[[i]] <- t_test
    i <- i + 1
  }
}

fc <- t_tests %>% lapply(., function(x) x$estimate["mean of x"] / x$estimate["mean of y"]) %>% unlist
pv <- t_tests %>% lapply(., function(x) x$p.val) %>% unlist
qv <- p.adjust(p = pv, method = "fdr")
mlog10Q <- -log10(qv)
mlog10P <- -log10(pv)

#######################################################################################
## Compute data for volcano plot
#######################################################################################

volcano_data <- data.frame(name = protein_names, log2FC = log2(unlist(fc)), mlog10Qval = mlog10Q %>% as.numeric, mlog10Pval = mlog10P,
                           comparison = rep(paste0(comparisons[1,], "vs", comparisons[2,]), each = length(protein_names))
)

fc_criteria <- abs(volcano_data$log2FC) > log.fold.change.threshold
pval_criteria <- ifelse(is.na(volcano_data$mlog10Pval > p.value.threshold), F, volcano_data$mlog10Pval > p.value.threshold)
both_criteria <- fc_criteria & pval_criteria
none_criteria <- !fc_criteria & !pval_criteria

volcano_data[fc_criteria, "criteria"] <- "FC"
volcano_data[pval_criteria, "criteria"] <- "PVAL"
volcano_data[both_criteria, "criteria"] <- "BOTH"
volcano_data[none_criteria, "criteria"] <- "NONE"
volcano_data$criteria <- volcano_data$criteria %>% factor(levels = c("FC", "PVAL", "BOTH","NONE"))
palette[which(levels(volcano_data$criteria) %in% dep_criteria)] <- "green"

p <- ggplot(volcano_data, aes(x = log2FC, y = mlog10Pval, col = criteria)) + geom_point() +
  geom_text(data = volcano_data %>% filter(criteria == "BOTH"),
           mapping = aes(x = log2FC, y = mlog10Pval, label = name)) +
  geom_hline(yintercept = p.value.threshold, linetype = "dashed") +
  geom_vline(xintercept = c(-log.fold.change.threshold, log.fold.change.threshold), linetype = "dashed") +
  labs(x = "log2(Fold change)", y = "-log10(P)") +
  scale_color_manual(values = palette) +
  guides(col=F) +
  ggtitle("Volcano plot")
p
ggsave(filename = file.path(output_dir, "volcano_plot.png"), plot = p)


#######################################################################################
## Differentially expressed proteins
#######################################################################################
up_prots <- volcano_data %>% filter(criteria %in% dep_criteria, log2FC > 0) %>%
  .$name %>% as.character()

down_prots <- volcano_data %>%
  filter(criteria %in% dep_criteria, log2FC < 0) %>%
  .$name %>% as.character()

stable_prots <- volcano_data %>%
  filter(criteria == "PVAL") %>%
  .$name %>% as.character()

dap_prots <- c(up_prots, down_prots)


# library("gProfileR")
# gsea_pos <- gprofiler(query = up_prots, organism = "hsapiens")
# gsea_neg <- gprofiler(query = down_prots, organism = "hsapiens")
# gsea_stable <- gprofiler(query = stable_prots, organism = "hsapiens")
# gsea <-  gprofiler(query = c(dap_prots, stable_prots), organism = "hsapiens")


write.table(x = volcano_data, file = file.path(output_dir, "volcano_data.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(x = up_prots, file = file.path(output_dir, "up_regulated.txt"), quote = F, row.names = F, col.names = F)
write.table(x = down_prots, file = file.path(output_dir, "down_regulated.txt"), quote = F, row.names = F, col.names = F)
write.table(x = dap_prots, file = file.path(output_dir, "dap_prots.txt"), quote = F, row.names = F, col.names = F)
write.table(x = stable_prots, file = file.path(output_dir, "stable_prots.txt"), quote = F, row.names = F, col.names = F)


#######################################################################################
## PCA
#######################################################################################

make_pca <- function(quant) {

  protein_names <- colnames(quant)
  colnames(quant) <- NULL
  
  prcomp_res <- prcomp(x = quant, scale. = T)
  prcomp_res_x <- prcomp_res$x %>% as.data.frame
  prcomp_res_x$condition <- substr(rownames(prcomp_res_x), 1,4)


  # prcomp_res_x <- full_join(prcomp_res_x, data %>% select(file, rep, condition) %>% distinct)
  
  p <- ggplot(data = prcomp_res_x, mapping = aes(
    col=factor(condition),
    x = PC1, y = PC2)
    ) +
    geom_point(size = 2) +
    #geom_text(aes(label = replicate))
    # scale_color_viridis(discrete = T) +
    # facet_grid(. ~ V1) +
    ggtitle(label = "PCA")
  ggsave(filename = file.path(output_dir, "pca.png"), plot = p)
  
  return(p)
}

quant <- data_spread %>% as.matrix %>% t
make_pca(quant)



# #######################################################################################
# ## Cumulative variance across PCs
# #######################################################################################
# 
# cumvar_data <- ((prcomp_res$sdev %>% cumsum) / sum(prcomp_res$sdev)) %>%
#   data.frame(PC=1:length(.), variance=.)
# 
# p <- ggplot(cumvar_data, mapping = aes(x = PC, y = variance)) +
#   geom_point()
# p
# ggsave2(filename = "cumulative_variance.png", plot = p)
# 
# 
# #######################################################################################
# ## Compute protein co expression network
# #######################################################################################
# 
# data_average <- data %>% na.omit %>% select(name, condition, Quantification) %>%
#   group_by(condition, name) %>%
#   summarise(average = mean(Quantification)) %>%
#   spread(name, average) %>%
#   column_to_rownames(var="condition") %>%
#   as.matrix
# data_average[is.na(data_average)] <- 0
# 
# cor_matrix <- cor(data_average)
# sign_matrix <- cor_matrix > 0
# 
# thresholds <- seq(from=0.80,to=0.99,by=0.01)
# mean.connectivities <- vector(length=length(thresholds))
# 
# 
# scale.free.R2 <- vector(length=length(thresholds))
# library(igraph)
# #3.3.-Espectro exhaustivo de posibles redes en funci?n de los umbrales
# #Itera para cada valor de umbral potencial de correlaci?n.
# for(i in 1:length(thresholds))
# {
#   print(thresholds[i])
#   ## Las siguientes instrucciones construyen la red basada en el umbral actual
#   ## Primero necesitamos crear la matriz de adyacencia, que es la representaci?n
#   ## matem?tica de cualquier red
#   
#   ## MATRIZ DE ADYACENCIA
#   ## La matriz de adyacencia se construye realizando la comparaci?n entre
#   ## cada elemento de la matriz de correlaci?n y el umbral de correlaci?n actual
#   ## de manera que aquellas posiciones ij de la matriz de correlacion que sean
#   ## mayores que el umbral dar?n un TRUE en la misma posicion de la matriz de adyacencia
#   ## y aquellas posiciones con un valor menor que el umbral, dar?n un FALSE. Como no quiero
#   ## conectar un gen consigo mismo, elimino los TRUEs producidos por la correlaci?n de un gen
#   ## consigo mismo, es decir, las correlaciones cuyo valor sea 1
#   current.adjacency <- (cor_matrix > thresholds[i] & cor_matrix < 1)
#   ## Para comprobar lo que el c?digo est? haciendo imprimiremos una posici?n aleatoria de la matriz de adyacencia actual
#   print(current.adjacency[3,4])
#   ## Tambi?n imprimimos la correlaci?n de un par de genes que ya sab?amos que estaban correlacionados,
#   ## con lo que esperamos un TRUE
#   print(current.adjacency[3,1])
#   ## We can see that the values stored in these matrix are logic boolean values: either TRUE or FALSE, depending
#   ## on whether the value stored in that same position of correlation was bigger or smaller than the current
#   ## threshold.
#   
#   ## RED A PARTIR DE LA MATRIZ DE ADYACENCIA
#   ## Build the network based on the TRUE and FALSE matrix (current.adjacency)
#   threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
#   
#   ## Compute node degrees. Se genera un vector num?rico en el que cada elemento corresponde
#   ## a un gen. Siempre es de la misma longitud, y esta longitud ser? igual a la cantidad de DEGS
#   ## obtenidos. Si un gen no tiene aristas no se elimina, sino que tendr? grado 0.
#   ## El vector est? nombrado, cada elemento tiene el nombre del gen cuyo grado muestra
#   node.degrees <- degree(threshold.network)
#   
#   ## Keep track of the mean connectivity or mean node degree of the current network
#   mean.connectivities[i] <- mean(node.degrees)
#   
#   ## Check scale free property
#   h <- hist(node.degrees)
#   ## El histograma se parecer? m?s a una exponencial negativa cuanto
#   ## m?s libre de escala sea la red lo que consigo subiendo el umbral de corte
#   
#   ## Compute degree frequencies, la frecuencia de cada uno de los grados.
#   ## La tabla recoge la frecuencia absoluta de cada grado, es decir, ser? una tabla en la
#   ## que cada elemento est? formado por dos n?meros: abajo se muestra la cantidad de nodos
#   ## que tienen el grado mostrado arriba.
#   degree.frequencies <- table(node.degrees)
#   
#   ## Determine linear regression for logarithmic transformed degree frequencies. Este paso
#   ## me permitir? extraer la r? asociada al ajuste de los datos de distribuci?n de grado de nodos
#   ##a una exponencial negativa
#   lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
#   ## Extract R squared as a measure of adjustment to the scale free property
#   s.lm <- summary(lm.r)
#   scale.free.R2[i] <- s.lm[["adj.r.squared"]]
# }
# 
# names(mean.connectivities) <- thresholds
# mean.connectivities
# ## El elemento i del vector scale.free.R2 guarda el coeficiente r? que mide la bondad del ajuste
# ## de la distribuci?n del grado de los nodos de la red i a una exponencial negativa.
# ## Cuanto m?s grande sea este par?metro, mi red ser? m?s cercana a una red libre de escala pura
# scale.free.R2
# 
# par(mfrow=c(1,2))
# plot(thresholds,mean.connectivities,type="o",col="red",
#      lwd=3,xlab="Correlation Threshold", ylab="Mean connectivity")
# plot(thresholds,scale.free.R2,type="o",col="blue",lwd=3, xlim=c(0.80,0.99),
#      xlab="Correlation Threshold", ylab="Scale Free Model Fit (R2)")
# par(mfrow=c(1,1))
#   
# adjacency_matrix <- cor_matrix > correlation.threshold
# 
# network <- graph.adjacency(adjacency_matrix, mode="undirected")
# 
# write.graph(network, file=file.path(output_dir, "network.gml"), format="gml")
# 
# 
# theme_networkMap <- theme(
#   #plot.background = element_rect(fill = "beige"),
#   panel.border = element_blank(),
#   panel.grid = element_blank(),
#   panel.background = element_blank(),
#   legend.background = element_blank(),
#   legend.position = "none",
#   legend.title = element_text(colour = "black"),
#   legend.text = element_text(colour = "black"),
#   legend.key = element_blank(),
#   axis.text = element_blank(), 
#   axis.title = element_blank(),
#   axis.ticks = element_blank()
# )
# 
# library("ggraph")
# library("tidygraph")
# 
# # VISUALIZE NETWORK
# p <- ggraph(network %>% as_tbl_graph %>%
#          mutate(Popularity = centrality_degree(mode = 'in')), layout = "auto") +
#   geom_edge_fan(aes(alpha = ..index..), show.legend = FALSE) + 
#   geom_node_point(aes(size = Popularity)) + 
#   # geom_edge_diagonal(label_colour = "blue") +
#   geom_node_point() +
#   theme_networkMap
# 
#   # scale_fill_gradient(high = "blue", low = "lightblue")
#   # labs(title = "Coauthorship Network of Jaap Paauwe",
#   #      subtitle = "Publications with more than one Google Scholar citation included",
#   #      caption = "paulvanderlaken.com") 
# ggsave2("network.png", plot = p)

