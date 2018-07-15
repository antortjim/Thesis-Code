library("ggplot2")
library("dplyr")
library("stringr")
library("viridis")
library("cowplot")
library(Cairo)
theme_set(theme_bw())
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
exp_dir <- "thesis/genedata/maxlfq"
input_dir <- file.path(home_dir, exp_dir, "/peptideShaker_out/reports")
#reports <- read.table(file.path(input_dir, "Default_PSM_Report.txt"), sep = "\t", header=T)
sample_name <- "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_black_Frac13_3"

reports <- read.table(file.path(input_dir, paste0(sample_name, "_Extended_PSM_Report", ".txt")),
                      sep = "\t", header=T, stringsAsFactors = F)




reports$Validation <- factor(reports$Validation, levels = rev(c("Confident", "Doubtful", "Not Validated")))

colnames(reports)[19] <- "mz_error_ppm"
reports <- filter(reports, Decoy == 0)

ggplot(data = reports, aes(x=Measured.Charge, fill=Validation)) +
  geom_histogram(stat="count")

cols <- rev(c("#6ec461",
              # "#fcd322",
              "gray90"))

reports$Validate <- ifelse(reports$Validation == "Not Validated", "Not Validated", "Validated")


ggplot(data = reports, aes(x = mz_error_ppm, fill = Validate)) +
  geom_histogram() +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 1400, 200)) +
  labs(x="Precursor m/z error (ppm)", y = "Number of PSMs") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.box = "horizontal", panel.border = element_rect(colour = "black", fill=NA),
        legend.justification = "center") +
  ggtitle("QC Plot")

ggsave(filename = file.path(home_dir, exp_dir, "figures", "qc.png"), width=5.74, height=7)
ggsave(filename = file.path("C:/Users/aoj/Desktop/Thesis-Report/", "plots", "qc.png"),
       width=5.74, height=7)

n_spectra <- system(paste0("grep -c 'BEGIN IONS' ",
              file.path(home_dir, exp_dir, "data", "mgf", sample_name), ".mgf"), intern = T) %>%
  as.integer


stats <- reports %>% group_by(Validation) %>% summarise(count=n())

setwd("C:/Users/aoj/Desktop/peptideShaker")

## Read PS report
#####################################

custom_report <- read.table(file.path(input_dir, paste0(sample_name, "_Custom_PSM", ".txt")),
                      sep = "\t", header=T, stringsAsFactors = F)


custom_report <- custom_report[, -1]
column_names <- colnames(custom_report)[-1]
custom_report <- custom_report[, -ncol(custom_report)]
colnames(custom_report) <- column_names

# custom_report %>% filter(Spectrum.Title == "mzspec:PXD000279:20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_black_Frac13_3.RAW:scan:10000")
                         
custom_report[, "Algorithm.Score"] <- custom_report[, "Algorithm.Score"] %>%
  str_match(string = ., pattern = "MS-GF\\+ \\((\\d\\.\\d*.*)\\)") %>% .[,2] %>% as.numeric
colnames(custom_report)[grep(pattern = "Confidence", x = colnames(custom_report))] <- "Confidence"
colnames(custom_report)[grep(pattern = "Score", x = colnames(custom_report))] <- "Score"

custom_report$Decoy <- factor(ifelse(custom_report$Decoy == 0, "Target", "Decoy"),
                              levels=c("Target", "Decoy"))

custom_report <- select(custom_report, Protein.s., Measured.Charge, Spectrum.Title, RT, m.z, Score, Confidence, Decoy, Validation)

custom_report$mlogScore <- -log(custom_report$Score)

custom_report_best <- custom_report %>% group_by(Spectrum.Title, Decoy) %>%
  slice(which.max(mlogScore))

temp <- (custom_report_best$Spectrum.Title %>% table)
# random_duplicate <- sample(temp[temp > 2] %>% names, 1)
# custom_report_best %>% filter(Spectrum.Title == random_duplicate)



facet_names <- c(
  `TRUE` = "PSM >2+",
  `FALSE` = "PSM 2+"
)

ggplot(data = custom_report_best, aes(x = mlogScore, fill=factor(Decoy))) +
  geom_histogram(position="dodge", bins=20) +
  guides(fill=guide_legend(title="Database")) +
  scale_fill_manual(values = c("#7dcb72", "#ff0000")) + 
  facet_wrap(~(Measured.Charge != "2+"), labeller = as_labeller(facet_names)) +
  labs(y = "Frequency (absolute)")

#######

predictions <- custom_report %>% select(Score, Decoy)
predictions$pred <- predictions$Score * 1/max(predictions$Score)

source("https://raw.githubusercontent.com/joyofdata/joyofdata-articles/master/roc-auc/plot_pred_type_distribution.R")

predictions$survived <- ifelse(predictions$Decoy == "Target",0, 1)
plot_pred_type_distribution(predictions, 0.7)


## Find the score of the first validated hit.
## This value is the threshold set by PeptideShaker
custom_report_best <- custom_report_best %>% arrange(-Score)
xintercept <- custom_report[which(custom_report_best$Validation != "Not Validated")[1],]
xintercept <- xintercept$mlogScore

ggplot(data = custom_report_best %>% filter(Measured.Charge == "2+"),
       mapping = aes(x=mlogScore, y=Confidence)) + geom_line(col="blue") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by=5)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by=10)) +
  geom_vline(xintercept = -log(xintercept), color="red", size=1) +
  labs(x="Score", y = "Confidence  [%]") +
  ggtitle("Score vs. Confidence") +
  theme_bw()

ggsave("plot1_R.png")

### EVALUATE MATCHING PERFORMANCE
exp_design <- read.table(file = file.path(home_dir, exp_dir, "data", "experimental_design.tsv"),header=T)
count_spectra <- read.table(file = file.path(home_dir, exp_dir, "count_spectra.txt"), stringsAsFactors = F, col.names = c("Name", "total", "match"))
count_spectra$Name <- count_spectra$Name %>% strsplit(., split = "\\/") %>% lapply(., function(x) unlist(strsplit(rev(x)[[1]], split = "\\."))[1]) %>% unlist
count_spectra <- left_join(count_spectra, select(exp_design, Fraction, Experiment, Replicate, Name), by = "Name")

count_spectra$`% Matched` <- count_spectra$match / count_spectra$total

match_percent_rect <- ggplot(count_spectra, aes(xmin = (Fraction-.5), xmax = (Fraction+.5),
                                                ymin = Replicate-.5, ymax=Replicate+.5,
                          fill=`% Matched`)) +
  geom_rect() +
  facet_wrap(facets = ~Experiment, nrow = 2) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1,2,3)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_distiller(type="seq", palette = "Purples", direction = 1) +
  labs(x = "Fraction", y = "Replicate") +
  # coord_equal(ratio=1) +
  # scale_x_reverse(expand = c(0,0))
  theme(legend.position = "top", legend.box="horizontal", legend.justification = "center")

match_percent_rect

plot_legend <- get_legend(match_percent_rect)

match_percent_rect <- match_percent_rect + guides(fill=F)
match_percent_rect

ggsave(filename = file.path(home_dir, exp_dir, "figures", "match_percent_rect.png"))

count_spectra %>% summarise(mean = mean(`% Matched`))
count_spectra %>% filter(Fraction < 14) %>% summarise(mean = mean(`% Matched`))
count_spectra %>% filter(Fraction >= 14) %>% summarise(mean = mean(`% Matched`))

dat <- count_spectra %>% group_by(Fraction) %>% summarise(mean = mean(`% Matched`), sd = sd(`% Matched`))

ggplot(count_spectra, aes(y=`% Matched`,x=Fraction)) +
  geom_point(aes(col=factor(Replicate))) + geom_smooth()

match_percent_overplot <- ggplot(dat, aes(y=mean,x=Fraction)) +
  geom_point() + geom_line(col="blue") +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), alpha=.7, width=.5) +
  labs(y = "Mean matched [%]") +
  scale_x_continuous(expand = c(0, 0))
  
ggsave(plot = match_percent_overplot, filename = file.path(home_dir, exp_dir, "figures", "match_percent_overplot.png"))
match_percent_overplot

plot_list <- list(match_percent_rect, match_percent_overplot) %>%
  lapply(function(x) x #+ theme(plot.margin = unit(c(0,0,0,2), "lines"))
         )

combination <- plot_grid(plotlist = plot_list,
                         #rel_heights = c(0.3,0.7), rel_widths = c(1,1),
                         ncol=1, nrow=2, labels = c("A", "B"), vjust = c(1,1))
combination <- plot_grid(plot_legend,
          combination,
          nrow=2, rel_heights = c(0.1, 1))
ggsave(plot = combination, filename = file.path(home_dir, exp_dir, "figures", "match_percent.png"))
