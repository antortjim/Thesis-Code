library("ggplot2")
library("dplyr")
library("stringr")
library("viridis")
library("cowplot")
library("tidyr")
library(Cairo)
theme_set(theme_bw())
rm(list = ls())
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
exp_dir <- "thesis/genedata/maxlfq"
input_dir <- file.path(home_dir, exp_dir, "/peptideShaker_out/PSM_reports/output_moff_RAW/mbr_output")
setwd(input_dir)

column_names <- c("run1", "run2", "run3")
combns <- combn(column_names, m = 2)
combns_names <- apply(combns, 2, function(x) paste(x, collapse = "vs"))

mbr_report_3 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_black_Frac13_3_match.txt", sep = "\t",header=T,stringsAsFactors = F)
mbr_report_2 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_green_Frac13_2_match.txt", sep = "\t",header=T,stringsAsFactors = F)
mbr_report_1 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_red_Frac13_1_match.txt", sep = "\t",header=T,stringsAsFactors = F)
mbr_report_1$replicate <- "run1"
mbr_report_2$replicate <- "run2"
mbr_report_3$replicate <- "run3"

mbr_report <- rbind(mbr_report_1, mbr_report_2, mbr_report_3)

mbr_summary <- group_by(mbr_report, replicate, matched) %>% summarise(count = n())

ggplot(mbr_summary, aes(x =replicate, y = count, fill=factor(matched))) + geom_bar(stat="identity", position="dodge")

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
  ifelse(., "MBR", "Identified")

colnames(matched_pairwise) <- combns_names
matched <- cbind(matched$peptide, matched_pairwise) %>% as.data.frame()

df <- df %>% select(-matched)
df <- df %>% spread(key = replicate, value = rt)

# df <- df[complete.cases(df),]
df <- cbind(df, matched)

diffs <- apply(combns, 2, function(x) df[,x[2]] - df[,x[1]]) %>% as.data.frame
colnames(diffs) <- combns_names

combns <- list(1:2, 3:4, 5:6) %>% lapply(function(x) combns[c(x[1], x[2])])


i <- 1
df_long <- data.frame(peptide=NULL, run=NULL, diff=NULL, Status=NULL, pair=NULL)

for(cmb in combns) {
  print(cmb)
  df_sub <- cbind(df[, c("peptide", cmb[1])], diffs[,i], matched_pairwise[, i], colnames(diffs)[i])
  colnames(df_sub)[2:5] <- c("run", "diff", "Status", "pair")
  df_long <- rbind(df_long, df_sub)
  i <- i + 1
}


p <- ggplot(df_long %>% filter(abs(diff) < 3.5), aes(x=run, y=diff)) +facet_wrap(facets = ~pair, ncol=3)

if(color) {
  p <- p + geom_point(aes(col=Status), size=1)
} else {
  p <- p + geom_point()
}

p <- p + geom_hline(yintercept = 0, linetype="dashed") +
  viridis::scale_color_viridis(discrete=T)
p <- p + labs(x = "Retention time in first MS run [min]")
p <- p + labs(y = "Retention time diff between runs [min]")
p <- p + theme(legend.position = "top")
p

ggsave(filename = file.path("C:/Users/aoj/Desktop/Thesis-Report/", "plots", "mbr.png"),
       width=5.74, height=7)


rm(list=ls())
home_dir <- ifelse(Sys.info()["sysname"] == "Windows", "//hest/aoj", "/z/home/aoj")
exp_dir <- "thesis/genedata/maxlfq"
input_dir <- file.path(home_dir, exp_dir, "/peptideShaker_out/PSM_reports/output_moff_RAW/mbr_output")

setwd(file.path(home_dir, exp_dir, "/peptideShaker_out/PSM_reports/output_moff_RAW/"))

apex_report_3 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_black_Frac13_3_match_moff_result.txt", sep = "\t",header=T,stringsAsFactors = F)
apex_report_2 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_green_Frac13_2_match_moff_result.txt", sep = "\t",header=T,stringsAsFactors = F)
apex_report_1 <- read.table(file = "20070904_CL_Orbi4_Offgel_XIC_Hela60_Ecoli10_red_Frac13_1_match_moff_result.txt", sep = "\t",header=T,stringsAsFactors = F)
apex_report_1$replicate <- "run1"
apex_report_2$replicate <- "run2"
apex_report_3$replicate <- "run3"
