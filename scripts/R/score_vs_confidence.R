library("ggplot2")
library("dplyr")
library("stringr")

## Read PS report
#####################################
filename <- "Custom_PSM.txt"
custom_report <- read.table(
  filename,
  stringsAsFactors = F, sep = "\t", fill=T, comment.char = "#", skip=1
)

column_names <- strsplit(readLines(filename,n = 1), split = "\\t") %>% unlist
column_names <- column_names %>% gsub(pattern = "\\#", replacement = "")
colnames(custom_report) <- column_names

custom_report[, "Algorithm Score"] <- custom_report[, "Algorithm Score"] %>%
  str_match(string = ., pattern = "MS-GF\\+ \\((\\d\\.\\d*.*)\\)") %>% .[,2] %>% as.numeric
custom_report <- custom_report %>% rename(score = `Algorithm Score`, confidence = `Algorithm Confidence [%]`)
custom_report <- custom_report[, -(c(1, ncol(custom_report)))]

## Process MS-GF+ scores
#####################################
custom_report$score <- -log(custom_report$score)
custom_report$score <- custom_report$score + 100 - max(custom_report$score)

## Find the score of the first validated hit.
## This value is the threshold set by PeptideShaker
custom_report <- custom_report %>% arrange(score)
xintercept <- custom_report[which(custom_report$Validation != "Not Validated")[1],"score"]

ggplot(data = custom_report %>% filter(`Identification Charge` == "2+"), aes(x=score, y=confidence)) + geom_line(col="blue") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by=5)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by=10)) +
  geom_vline(xintercept = xintercept, color="red", size=1) +
  labs(x="Score", y = "Confidence  [%]") +
  ggtitle("Score vs. Confidence") +
  theme_bw()

ggsave("plot1_R.png")