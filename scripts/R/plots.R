source("settings.R")

#######################################################################################
## Read data and save available conditions
#######################################################################################
data <- read.table(file = file.path(data_dir, "dataset.tsv"), stringsAsFactors = F, header = T, row.names = NULL)
data %>% group_by(condition1, condition2, rep) %>% summarise(count = n())

condition_fields <- grep(pattern = "condition", colnames(data), value = T)
conditions <- data[, condition_fields]

if(!"condition" %in% colnames(data)) {
data$condition <- apply(data[, condition_fields], 1, function(x) paste(x, collapse = "_"))
}

for (i in 1:length(analysis_conditions)) {
  data <- data[data[, names(analysis_conditions)[i]] %in% analysis_conditions[[names(analysis_conditions)[i]]],]
}
kept_files <- data$file %>% unique

#######################################################################################
## Historam of the counts of protein id occurrrences in replicates across conditions
#######################################################################################


counts_condition <- data %>%
  group_by(condition, name) %>% summarise(count = sum(!is.na(amount)))
p <- ggplot(counts_condition ,
       aes(x = count, fill=condition)) + geom_histogram(position="dodge") +
  scale_x_continuous(breaks = 0:max(counts_condition$count)) +
  labs(x="Available in # replicates", y = "Protein count")
ggsave2(filename = "counts_condition.png", plot = p)

# data <- data %>% filter(name %in% (data %>% group_by(condition, name) %>%
#                                      summarise(count = n()) %>%
#                                      filter(count > 4) %>% .$name))

# data[is.na(data$amount),"amount"] <- 0

#######################################################################################
## Boxplots of the protein amount distributions across samples
#######################################################################################

p <- ggplot(data %>% filter(!(name %in% excluded_proteins)),
       aes(x = rep, y = log(amount), fill=condition, group=rep)) +
  geom_boxplot() +
  facet_wrap(~condition) + guides(fill = F)
p
ggsave2(filename = "amount_boxplots.png", plot = p)

p <- ggplot(data %>% filter(!(name %in% excluded_proteins)), aes(x = condition, y = log(amount))) +
  geom_boxplot()
ggsave2(filename = "amounts_condition_boxplot.png", plot = p)

p <- ggplot(data %>% filter(!(name %in% excluded_proteins)), aes(x = amount, col = condition)) +
  geom_density(alpha = 0.5) + guides(col=F)
ggsave2(filename = "amount_density.png", plot = p)



#######################################################################################
## Scatterplot of the proteins amount across conditions (one datapoint per replicate)
#######################################################################################

ggplot(data %>% filter(!(name %in% excluded_proteins)), aes(x = name, y = amount, col=condition)) +
  geom_point() +
  facet_wrap(~condition) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#######################################################################################
## Scatter plot of mean and coefficient of variation
#######################################################################################

cv_data <- data %>%
  group_by(condition, name) %>%
  summarise(cv = sd(amount, na.rm = T) / mean(amount, na.rm = T),
            mu = mean(amount, na.rm=T),
            sigma = sd(amount, na.rm=T))

ggplot(cv_data, aes(x = mu, y = cv)) + geom_point() + facet_wrap(~condition)

# ggplot(cv_data %>% arrange(condition, -mu), aes(x = name, y = mu)) + geom_bar(stat="identity") + geom_errorbar(aes(x = name,
#                                                                     ymin = mu - sigma,
#                                                                     ymax = mu + sigma)) +
#   facet_wrap(~condition)
# 
# data$name %>% as.character %>% unique

#######################################################################################
## Proteins quantified per sample
#######################################################################################

quant_count <- data %>% na.omit %>% group_by(file) %>%
  summarise(count = n()) %>% mutate(file = as.factor(file))
ggplot(data = quant_count, aes(x = file, y = count)) + geom_bar(stat = "identity") + labs(y = "Protein count")


# data %>% filter(file == 1) %>% .$amount
# 
# join_data <- inner_join(data[data$file == 1, c("name", "amount")], data[data$file == 2, c("name", "amount")], by = "name")
# join_data[join_data %>% is.na] <- 0.1
# 
# M <- log2(join_data$amount.x) - log2(join_data$amount.y)
# A <- 0.5*(log2(join_data$amount.x) + log2(join_data$amount.y))
# 
# plot(A, M, pch = 16)


#######################################################################################
## Protein quantities across conditions
#######################################################################################

average_dat <- data %>%
  group_by(condition, name) %>% summarise(average = mean(amount)) %>% arrange(name)

average_dat_spread <- average_dat %>% spread(condition, average)
# average_dat_spread[complete.cases(average_dat_spread),]

p <- ggplot(average_dat, aes(x = condition, y = log(average), group = name)) +
  geom_line(alpha = 0.2) +
  guides(group = F, col = F)
ggsave2(filename = "logAverage_condition.png", plot = p)

# progression <- data_average %>% as.data.frame %>%
#   rownames_to_column(var = "condition") %>%
#   gather(name, amount, -condition)

protein_name <- sample(average_dat$name, 5)

ggplot(average_dat %>% filter(name %in% protein_name),
       aes(x = condition, y = average, col=name, group=name)) +
  geom_line() +
  guides(fill=F)



## Remove NAs
# unique_proteins_conditions <- analysis_conditions %>% lapply(function(x) {
#     average_dat_spread[(average_dat_spread %>% is.na %>% .[,x]),]$name
#   }
# )

average_dat_spread <- average_dat_spread[complete.cases(average_dat_spread),]


## Keep proteins that are quantified in all samples
protein_names <- average_dat_spread$name

## Compute possible comparisons
comparisons <- combn(unique(data$condition), m = 2)

## Compute fold change and significance
fc <- 1:ncol(comparisons) %>% lapply(function(i) {
  average_dat_spread[, comparisons[1, i]] / average_dat_spread[, comparisons[2, i]]
  })
p_values <- matrix(ncol = ncol(comparisons), nrow = length(protein_names))

for (j in 1:ncol(comparisons)) {
  i <- 1
  for(protein in protein_names) {
    #print(protein)
    commercial <- data %>% filter(condition == comparisons[1, j], name == protein) %>% .$amount
    purified <- data %>% filter(condition == comparisons[2, j], name == protein) %>% .$amount
    p_value <- tryCatch({
      t.test(commercial, purified, alternative = "two.sided") %>% .$p.value},
      error = function(e) {
        NA
    })
    p_values[i, j] <- p_value
    i <- i + 1
  }
}
mlog10P <- -log10(p_values)

#######################################################################################
## Compute data for volcano plot
#######################################################################################

volcano_data <- data.frame(name = protein_names, log2FC = log2(unlist(fc)), mlog10Pval = mlog10P %>% as.numeric,
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

p <- ggplot(volcano_data %>% filter(comparison == compare), aes(x = log2FC, y = mlog10Pval, col = criteria)) + geom_point() +
  geom_text(data = volcano_data %>% filter(criteria == "BOTH"),
           mapping = aes(x = log2FC, y = mlog10Pval, label = name)) +
  geom_hline(yintercept = p.value.threshold, linetype = "dashed") +
  geom_vline(xintercept = c(-log.fold.change.threshold, log.fold.change.threshold), linetype = "dashed") +
  labs(x = "log2(Fold change)", y = "-log10(P)") +
  scale_color_manual(values = palette) +
  guides(col=F) +
  ggtitle("Volcano plot")
p
ggsave2(filename = "volcano_plot.png", plot = p)


#######################################################################################
## Differentially expressed proteins
#######################################################################################
pos_deps <- volcano_data %>% filter(criteria %in% dep_criteria, comparison == compare, log2FC > 0) %>%
  .$name %>% as.character()

neg_deps <- volcano_data %>%
  filter(criteria %in% dep_criteria, comparison == compare, log2FC < 0) %>%
  .$name %>% as.character()

deps <- c(pos_deps, neg_deps)


library("gProfileR")

gsea <- gprofiler(query = deps, organism = "hsapiens")

#######################################################################################
## PCA
#######################################################################################

quant <- data %>%
  # keep the protein name, its quantity and the file (sample id)
  select(name, amount, file) %>%
  # make the data wide
  spread(key = name, value = amount) %>%
  select(-file) %>%
  as.matrix
colSums(quant %>% is.na) %>% table
quant <- quant[,colSums(quant %>% is.na) < (nrow(quant) / 5)]

## Imputation of missing values

col_means <- colMeans(quant, na.rm = T)
quant2 <- 1:ncol(quant) %>% lapply(function(i) ifelse(is.na(quant[,i]), col_means[i], quant[,i])) %>% do.call(cbind, .)
colnames(quant2) <- colnames(quant)
quant <- quant2


make_pca <- function(quant) {

  protein_names <- colnames(quant)
  colnames(quant) <- NULL
  
  prcomp_res <- prcomp(x = quant)
  prcomp_res_x <- prcomp_res$x %>% as.data.frame
  prcomp_res_x$file <- kept_files
  
  prcomp_res_x <- full_join(prcomp_res_x, data %>% select(file, rep, condition) %>% distinct)
  
  p <- ggplot(data = prcomp_res_x, mapping = aes(
    col=factor(condition),
    x = PC1, y = PC2)
    ) +
    geom_point(size = 2) +
    #geom_text(aes(label = replicate))
    scale_color_viridis(discrete = T) +
    # facet_grid(. ~ V1) +
    ggtitle(label = project_name %>% simpleCap)
  ggsave2(filename = "pca.png", plot = p)
  
  return(p)
}

make_pca(quant)



#######################################################################################
## Cumulative variance across PCs
#######################################################################################

cumvar_data <- ((prcomp_res$sdev %>% cumsum) / sum(prcomp_res$sdev)) %>%
  data.frame(PC=1:length(.), variance=.)

p <- ggplot(cumvar_data, mapping = aes(x = PC, y = variance)) +
  geom_point()
p
ggsave2(filename = "cumulative_variance.png", plot = p)


#######################################################################################
## Compute protein co expression network
#######################################################################################

data_average <- data %>% na.omit %>% select(name, condition, amount) %>%
  group_by(condition, name) %>%
  summarise(average = mean(amount)) %>%
  spread(name, average) %>%
  column_to_rownames(var="condition") %>%
  as.matrix
data_average[is.na(data_average)] <- 0

cor_matrix <- cor(data_average)
sign_matrix <- cor_matrix > 0

thresholds <- seq(from=0.80,to=0.99,by=0.01)
mean.connectivities <- vector(length=length(thresholds))


scale.free.R2 <- vector(length=length(thresholds))
library(igraph)
#3.3.-Espectro exhaustivo de posibles redes en funci?n de los umbrales
#Itera para cada valor de umbral potencial de correlaci?n.
for(i in 1:length(thresholds))
{
  print(thresholds[i])
  ## Las siguientes instrucciones construyen la red basada en el umbral actual
  ## Primero necesitamos crear la matriz de adyacencia, que es la representaci?n
  ## matem?tica de cualquier red
  
  ## MATRIZ DE ADYACENCIA
  ## La matriz de adyacencia se construye realizando la comparaci?n entre
  ## cada elemento de la matriz de correlaci?n y el umbral de correlaci?n actual
  ## de manera que aquellas posiciones ij de la matriz de correlacion que sean
  ## mayores que el umbral dar?n un TRUE en la misma posicion de la matriz de adyacencia
  ## y aquellas posiciones con un valor menor que el umbral, dar?n un FALSE. Como no quiero
  ## conectar un gen consigo mismo, elimino los TRUEs producidos por la correlaci?n de un gen
  ## consigo mismo, es decir, las correlaciones cuyo valor sea 1
  current.adjacency <- (cor_matrix > thresholds[i] & cor_matrix < 1)
  ## Para comprobar lo que el c?digo est? haciendo imprimiremos una posici?n aleatoria de la matriz de adyacencia actual
  print(current.adjacency[3,4])
  ## Tambi?n imprimimos la correlaci?n de un par de genes que ya sab?amos que estaban correlacionados,
  ## con lo que esperamos un TRUE
  print(current.adjacency[3,1])
  ## We can see that the values stored in these matrix are logic boolean values: either TRUE or FALSE, depending
  ## on whether the value stored in that same position of correlation was bigger or smaller than the current
  ## threshold.
  
  ## RED A PARTIR DE LA MATRIZ DE ADYACENCIA
  ## Build the network based on the TRUE and FALSE matrix (current.adjacency)
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  ## Compute node degrees. Se genera un vector num?rico en el que cada elemento corresponde
  ## a un gen. Siempre es de la misma longitud, y esta longitud ser? igual a la cantidad de DEGS
  ## obtenidos. Si un gen no tiene aristas no se elimina, sino que tendr? grado 0.
  ## El vector est? nombrado, cada elemento tiene el nombre del gen cuyo grado muestra
  node.degrees <- degree(threshold.network)
  
  ## Keep track of the mean connectivity or mean node degree of the current network
  mean.connectivities[i] <- mean(node.degrees)
  
  ## Check scale free property
  h <- hist(node.degrees)
  ## El histograma se parecer? m?s a una exponencial negativa cuanto
  ## m?s libre de escala sea la red lo que consigo subiendo el umbral de corte
  
  ## Compute degree frequencies, la frecuencia de cada uno de los grados.
  ## La tabla recoge la frecuencia absoluta de cada grado, es decir, ser? una tabla en la
  ## que cada elemento est? formado por dos n?meros: abajo se muestra la cantidad de nodos
  ## que tienen el grado mostrado arriba.
  degree.frequencies <- table(node.degrees)
  
  ## Determine linear regression for logarithmic transformed degree frequencies. Este paso
  ## me permitir? extraer la r? asociada al ajuste de los datos de distribuci?n de grado de nodos
  ##a una exponencial negativa
  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  ## Extract R squared as a measure of adjustment to the scale free property
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

names(mean.connectivities) <- thresholds
mean.connectivities
## El elemento i del vector scale.free.R2 guarda el coeficiente r? que mide la bondad del ajuste
## de la distribuci?n del grado de los nodos de la red i a una exponencial negativa.
## Cuanto m?s grande sea este par?metro, mi red ser? m?s cercana a una red libre de escala pura
scale.free.R2

par(mfrow=c(1,2))
plot(thresholds,mean.connectivities,type="o",col="red",
     lwd=3,xlab="Correlation Threshold", ylab="Mean connectivity")
plot(thresholds,scale.free.R2,type="o",col="blue",lwd=3, xlim=c(0.80,0.99),
     xlab="Correlation Threshold", ylab="Scale Free Model Fit (R2)")
par(mfrow=c(1,1))
  
adjacency_matrix <- cor_matrix > correlation.threshold

network <- graph.adjacency(adjacency_matrix, mode="undirected")

write.graph(network, file=file.path(output_dir, "network.gml"), format="gml")


theme_networkMap <- theme(
  #plot.background = element_rect(fill = "beige"),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  legend.background = element_blank(),
  legend.position = "none",
  legend.title = element_text(colour = "black"),
  legend.text = element_text(colour = "black"),
  legend.key = element_blank(),
  axis.text = element_blank(), 
  axis.title = element_blank(),
  axis.ticks = element_blank()
)

library("ggraph")
library("tidygraph")

# VISUALIZE NETWORK
p <- ggraph(network %>% as_tbl_graph %>%
         mutate(Popularity = centrality_degree(mode = 'in')), layout = "auto") +
  geom_edge_fan(aes(alpha = ..index..), show.legend = FALSE) + 
  geom_node_point(aes(size = Popularity)) + 
  # geom_edge_diagonal(label_colour = "blue") +
  geom_node_point() +
  theme_networkMap

  # scale_fill_gradient(high = "blue", low = "lightblue")
  # labs(title = "Coauthorship Network of Jaap Paauwe",
  #      subtitle = "Publications with more than one Google Scholar citation included",
  #      caption = "paulvanderlaken.com") 
ggsave2("network.png", plot = p)

