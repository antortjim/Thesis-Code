library("ggplot2")
library("dplyr")
library("stringr")
library("viridis")
library("cowplot")
library("tidyr")
library(Cairo)
library(latex2exp)


theme_set(theme_bw())
rm(list = ls())
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj",
                   ifelse(Sys.info()["user"] == "aoj", "/z/home/aoj",
                          "/home/antortjim/"))

data
exp_dir <- "thesis/genedata/maxlfq"
input_dir <- file.path(home_dir, exp_dir, "/peptideShaker_out/PSM_reports/output_moff_RAW/mbr_output")
thesis_report_dir <- "C:/Users/aoj/OneDrive - Novozymes A S/Thesis-Report/"
setwd(input_dir)

column_names <- paste0("Run ", 1:3)

combns <- combn(column_names, m = 2)
combns_names <- apply(combns, 2, function(x) paste(x, collapse = " - "))

color <- T
mbr_report_3 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_black_Frac13_3_match.txt", sep = "\t",header=T,stringsAsFactors = F)
mbr_report_2 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_green_Frac13_2_match.txt", sep = "\t",header=T,stringsAsFactors = F)
mbr_report_1 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_red_Frac13_1_match.txt", sep = "\t",header=T,stringsAsFactors = F)
mbr_report_1$replicate <- column_names[1]
mbr_report_2$replicate <- column_names[2]
mbr_report_3$replicate <- column_names[3]

mbr_report <- rbind(mbr_report_1, mbr_report_2, mbr_report_3)

mbr_summary <- group_by(mbr_report, replicate, matched) %>% summarise(count = n())

sources <- c("MS1", "MBR")
mbr_summary$matched <- factor(ifelse(mbr_summary$matched == 0, sources[1], sources[2]),
                              levels = rev(sources))

run_labels <- 1:3
names(run_labels) <- column_names

ggplot(mbr_summary, aes(x = replicate, y = count/1e3, fill=matched)) +
  geom_bar(stat="identity") +
  scale_fill_viridis(discrete = T) +
  guides(fill=guide_legend(title="Source")) +
  labs(x="Run", y=TeX("# 10^3 Matches")) +
  scale_x_discrete(labels=run_labels) +
  theme_bw(base_size=20)



ggsave(filename = file.path(thesis_report_dir, "plots", "mbr_summary.png"),
       width=6, height=3)


df1 <- select(mbr_report_1, peptide, rt, replicate, matched) %>% arrange(rt)
df1 <- df1[!duplicated(df1$peptide),]
df2 <- select(mbr_report_2, peptide, rt, replicate,matched) %>% arrange(rt)
df2 <- df2[!duplicated(df2$peptide),]
df3 <- select(mbr_report_3, peptide, rt, replicate,matched) %>% arrange(rt)
df3 <- df3[!duplicated(df3$peptide),]


df <- rbind(df1, df2, df3) %>%
  # filter(peptide != "") %>%
  arrange(peptide, replicate) %>%
  mutate(rt = rt/60)

matched <- df %>% select(-rt) %>% spread(key = replicate, value = matched)

matched_pairwise <- apply(combns, 2, function(x) matched[, x[1]] | matched[, x[2]]) %>%
  ifelse(., sources[2], sources[1])

colnames(matched_pairwise) <- combns_names
matched <- cbind(matched$peptide, matched_pairwise) %>% as.data.frame()

df <- df %>% select(-matched)
df <- df %>% spread(key = replicate, value = rt)

# df <- df[complete.cases(df),]
df <- cbind(df, matched)

diffs <- apply(combns, 2, function(x) df[,x[1]] - df[,x[2]]) %>% as.data.frame
colnames(diffs) <- combns_names

combns <- list(1:2, 3:4, 5:6) %>% lapply(function(x) combns[c(x[1], x[2])])


i <- 1
df_long <- data.frame(peptide=NULL, run=NULL, diff=NULL, Source=NULL, pair=NULL)

for(cmb in combns) {
  print(cmb)
  df_sub <- cbind(df[, c("peptide", cmb[1])], diffs[,i], matched_pairwise[, i], colnames(diffs)[i])
  colnames(df_sub)[2:5] <- c("run", "diff", "Source", "pair")
  df_long <- rbind(df_long, df_sub)
  i <- i + 1
}



p <- ggplot(df_long %>% filter(abs(diff) < 3.5), aes(x=run, y=diff)) + 
  facet_wrap(
    facets = ~pair, ncol=3
  )
p

if(color) {
  p <- p + geom_point(aes(col=Source), size=.05)
} else {
  p <- p + geom_point()
}

p <- p + geom_hline(yintercept = 0, linetype="dashed") +
  viridis::scale_color_viridis(discrete=T)
p <- p + labs(x = "RT in first MS run [min]")
p <- p + labs(y = "RT diff btw. runs [min]")
p <- p + theme_bw(base_size=20)
p <- p + guides(col = guide_legend(override.aes = list(size = 5)))
p <- p + theme(legend.position = "top")
p <- p + coord_cartesian(ylim = c(-2,2))
p
ggsave(filename = file.path(thesis_report_dir, "plots", "mbr.png"),
       width=6, height=5)


rm(list=ls())
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
exp_dir <- "thesis/genedata/maxlfq"
input_dir <- file.path(home_dir, exp_dir, "/peptideShaker_out/PSM_reports/output_moff_RAW/mbr_output")
thesis_report_dir <- "C:/Users/aoj/OneDrive - Novozymes A S/Thesis-Report/"
source(file.path(home_dir, "thesis", "genedata", "scripts", "R", "check_organism.R"))


setwd(file.path(home_dir, exp_dir, "/peptideShaker_out/PSM_reports/output_moff_RAW/"))
# 
# apex_report_3 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_black_Frac13_3_match_moff_result.txt", sep = "\t",header=T,stringsAsFactors = F)
# apex_report_2 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_green_Frac13_2_match_moff_result.txt", sep = "\t",header=T,stringsAsFactors = F)
# apex_report_1 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_red_Frac13_1_match_moff_result.txt", sep = "\t",header=T,stringsAsFactors = F)
# apex_report_1$replicate <- "run1"
# apex_report_2$replicate <- "run2"
# apex_report_3$replicate <- "run3"

exp_design <- read.table(file = file.path(home_dir, exp_dir, "data", "experimental_design.tsv"),header=T)
peptide_summary <- read.table(file = "peptide_summary_intensity_moFF_run.tab", sep = "\t",header=T,stringsAsFactors = F)
peptide_summary_long <- peptide_summary %>% gather(key = "Name", value = "apex_intensity", -peptide, -prot)
peptide_summary_long$Name <- peptide_summary_long$Name %>% gsub(pattern = "sumIntensity_", replacement = "", x = .)

peptide_summary_long <- left_join(peptide_summary_long, select(exp_design, Name, Experiment, Replicate, Fraction), by = "Name") %>%
  select(-Name)

# df1 <-   group_by(peptide_summary_long, Experiment, Replicate, Fraction) %>% summarise(count = sum(apex_intensity==0))
# df1$Status <- "missing"

df <- group_by(peptide_summary_long, Experiment, Replicate, Fraction) %>% summarise(q = 100*sum(apex_intensity!=0)/length(apex_intensity))
df <- group_by(peptide_summary_long, Experiment, Replicate, Fraction) %>% summarise(q = sum(apex_intensity!=0))

# df2$Status <- "available"

# df <- rbind(df1, df2)

ggplot(data = df, aes(xmin=Fraction-0.5, xmax=Fraction+0.5,
                      ymin=Replicate-0.5, ymax=Replicate+0.5,
                      fill=q)) +
  facet_wrap(~Experiment, nrow=2) +
  geom_rect() +
  scale_fill_distiller(direction=1)


l_samples <- grep(pattern = "Ecoli10", colnames(peptide_summary))
h_samples <- grep(pattern = "Ecoli30", colnames(peptide_summary))


x <- peptide_summary[
  rowSums(peptide_summary[,h_samples] == 0) > 3 &
  rowSums(peptide_summary[,l_samples] == 0) > 3
    ,"peptide"]

ratio <- (peptide_summary[,h_samples] %>% apply(.,1,max))  / (peptide_summary[,l_samples] %>% apply(.,1,max))
x_human <- ratio > 0.5 & ratio < 2 & !is.infinite(ratio)
x_ecoli <- ratio > 2 & !is.infinite(ratio)


peptide_human <- sample(peptide_summary[x_human, "peptide"], 1)
prot_human <- filter(peptide_summary, peptide == peptide_human) %>% .$prot %>%
  strsplit(., split = ", ") %>% unlist %>% .[1]

p <- ggplot(data = peptide_summary_long %>% filter(peptide == peptide_human),
       aes(Fraction, y = apex_intensity, group=Replicate)) +
  facet_wrap(~Experiment, nrow=2) +
  geom_line(col="blue") +
  labs("MS1 Apex intensity") +
  ggtitle(peptide_human, subtitle = paste0(prot_human, ", Homo sapiens"))

peptide_ecoli <- sample(peptide_summary[x_ecoli, "peptide"], 1)
prot_ecoli <- filter(peptide_summary, peptide == peptide_ecoli) %>% .$prot %>%
  strsplit(., split = ", ") %>% unlist %>% .[1]

q <- ggplot(data = peptide_summary_long %>% filter(peptide == peptide_ecoli),
       aes(Fraction, y = apex_intensity, group=Replicate)) +
  facet_wrap(~Experiment, nrow=2) +
  geom_line(col="red") +
  ggtitle(peptide_ecoli, subtitle = paste0(prot_ecoli, ", E. coli"))

pp <- plot_grid(p+labs(y="Apex intensity"),q+labs(y=""))
pp
ggsave(plot = pp, filename = file.path(thesis_report_dir, "plots", "peptide_profile.png"),
       width=7, height=4)

human_maxI <- cbind(peptide_summary[x_human,1:2], ratio_apex_intensity = 
                    (peptide_summary[x_human,h_samples] %>% apply(.,1,max)) / (peptide_summary[x_human,l_samples] %>% apply(.,1,max))
)

ecoli_maxI <- cbind(peptide_summary[x_ecoli,1:2], ratio_apex_intensity = 
                      (peptide_summary[x_ecoli,h_samples] %>% apply(.,1,max)) / (peptide_summary[x_ecoli,l_samples] %>% apply(.,1,max))
                    
      )

df <- rbind(cbind(human_maxI, taxon="Homo sapiens"),
      cbind(ecoli_maxI, taxon="E. coli")
      )

df$taxon <- factor(as.character(df$taxon), levels = c("E. coli", "Homo sapiens"))

ratio_boxplot <- boxplot(df$ratio_apex_intensity)
outlier_min <- ratio_boxplot$out %>% min

pp <- ggplot(data = df %>% filter(ratio_apex_intensity < outlier_min),
       aes(x=ratio_apex_intensity, fill=taxon)) +
  geom_density(alpha=0.5) +
  labs(x = "Apex intensity ratio [H/L]", y = "Density") +
  guides(fill=guide_legend(title="Organism")) + theme(legend.position="bottom") +
  scale_x_continuous(breaks = seq(0, ceiling(outlier_min), 1))
pp

ggsave(plot = pp, filename = file.path(thesis_report_dir, "plots", "density_ratio.png"),
       width=7, height=4)


peptide_summary_long$prot %>% unique %>% length
organism_check <- check_organism(proteins = peptide_summary_long$prot %>% unique,split = ", ")

taxonomy <- data.frame(taxon = organism_check[[2]], prot = unique(peptide_summary_long$prot))
peptide_summary_long <- left_join(peptide_summary_long, taxonomy, by="prot")

peptide_summary_long <- peptide_summary_long %>% filter(!is.na(taxon))

filter_peps <- peptide_summary_long %>% group_by(Experiment, peptide) %>% summarise(count=n())
filter_peps <- filter_peps %>% select(-count)
peptide_h <- filter_peps %>% filter(Experiment == "H") %>% .$peptide
peptide_l <- filter_peps %>% filter(Experiment == "L") %>% .$peptide

peptides <- peptide_h[peptide_h %in% peptide_l]


plot_data <- peptide_summary_long %>% filter(peptide %in% peptides)

plot_data %>% group_by(Experiment, taxon) %>% summarise(count = n()) %>% arrange(taxon)

ggplot(plot_data, aes(x=apex_intensity, fill = Experiment)) +
  # geom_histogram(bins=100)
  geom_density(alpha=0.5) +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~taxon)
  

