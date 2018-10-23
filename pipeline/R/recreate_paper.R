library("cowplot")
library("ggplot2")
theme_set(theme_bw())
library("optparse")

home_dir <- ifelse(Sys.info()["sysname"] == "Linux", "/z/home/aoj", "//hest/aoj")
input_dir <- file.path(home_dir, "/thesis/genedata/maxlfq/paper")
option_list = list(
  make_option(c("--root_dir"), type="character", default=file.path(home_dir, "/thesis/genedata")),
  make_option(c("--exp_name"), type="character", default="maxlfq"),
  make_option(c("--input_dir"), type="character"),
  make_option(c("--output_dir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root_dir <- opt$root_dir
exp_name <- opt$exp_name


# Peptide intensities saved in proteinGroups.txt (not LFQ, raw intensities)
intensities <- read.table(file = file.path(input_dir, "intensities.tsv"), sep = "\t", header=T)
Protein.IDs <- read.table(file = file.path(input_dir, "Protein.IDs.tsv"), sep = "\t", header=T)[,1]
log2ratioLFQ <- read.table(file = file.path(input_dir, "log2ratioLFQ.tsv"), sep = "\t", header=T)[,1]
log2ratioIntensities <- read.table(file = file.path(input_dir, "log2ratioIntensities.tsv"), sep = "\t", header=T)[,1]
log2ratioSC <- read.table(file = file.path(input_dir, "log2ratioSC.tsv"), sep = "\t", header=T)[,1]

######################################################################
## Compute the sum of intensities for every protein
######################################################################

  ## For every protein i
    ## For every peptide j
       ## Sum the intensity in the 1:1 samples (i.e H1, H2, and H3)
    ## Compute the median of the sums 1:j
  ## The median of the sums is the summed intensity of protein i

summed_peptide_intensities <- numeric(length(Protein.IDs))
names(summed_peptide_intensities) <- Protein.IDs
reference_samples <- grep(pattern = "H", x = colnames(intensities))


#any(rowSums(supplementary[,1:3] == 0) + rowSums(supplementary[,4:6] == 0) > 2)

############################################################################################
## Create a data.frame containing all the summed intensities and the fold changes
############################################################################################

q_methods <- c("SC", "Intensity", "LFQ")
plot_data <- data.frame(
  log10sumI = rep(log10sumI, times=3),
  log2ratio = c(log2ratioSC, log2ratioIntensities, log2ratioLFQ),
  protein = rep(Protein.IDs, times=3),
  taxonomy = factor(
    rep(taxonomy, times=3),
    levels = c("Homo sapiens", "Escherichia coli (strain K12)")),
  method = factor(
    x = rep(q_methods, each = length(Protein.IDs)),
    levels = q_methods)
  )

plot_data %>% group_by(taxonomy, method)%>% summarise(count = n())


############################################################################################
## Recreate the MaxLFQ paper plots
############################################################################################

xbreaks <- seq(
  round(min(plot_data$log2ratio[!is.infinite(plot_data$log2ratio)])),
  round(max(plot_data$log2ratio[!is.infinite(plot_data$log2ratio)])),
  1
)

p <- ggplot(
  data = plot_data,
  mapping = aes(x = log2ratio, y = log10sumI, color = taxonomy)
  )  +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point() + facet_wrap(~method) +
  scale_x_continuous(breaks = xbreaks,limits = c(-2, 4))


Mode <- function(x, nbins) {
   x <- hist(x, breaks = nbins)
   x <- x$breaks[which.max(x$counts)]
   return(x)
  }

nbins <- 100
x <- hist(plot_data %>%
            filter(taxonomy == "Escherichia coli (strain K12)", method == "SC") %>%
            .$log2ratio, breaks = nbins, plot = F)
x$breaks[which.max(x$counts)]


modes <- plot_data %>% group_by(taxonomy, method)%>% summarise(mode = Mode(log2ratio, nbins)) %>%
  arrange(method, taxonomy) %>% .$mode
vlines <- data.frame(x = modes, method = rep(q_methods, each = 2))
  
lapply(c(0, 2, 4), function(x) modes[2+x] - modes[1+x]) %>% unlist

# ggplot(data = plot_data, mapping = aes(x = log2ratio, fill = taxonomy)) +
#   # geom_vline(data = vlines, mapping = aes(xintercept = x), alpha = 0.5, linetype="dashed") +
#   geom_histogram(bins=nbins, position="identity") +
#   facet_wrap(taxonomy~method)

q <- ggplot(data = plot_data, mapping = aes(x = log2ratio, fill = taxonomy)) +
  geom_vline(data = vlines, mapping = aes(xintercept = x), alpha = 0.5, linetype="dashed") +
  geom_histogram(bins=nbins, position = "identity") +
  facet_wrap(~method) +
  scale_x_continuous(breaks = xbreaks,limits = c(-2, 4))

paper_figure <- cowplot::plot_grid(p, q, nrow = 2, labels = c("A", "B"))

plot_data %>% group_by(taxonomy, method) %>%
  summarise(stdev = sd(log2ratio[!is.infinite(log2ratio)]))
