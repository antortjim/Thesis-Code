library("dplyr")
setwd("/z/home/aoj/thesis/genedata/thp1/")
msgfplus_output <- read.table("searchgui_out/PD7508-GDTHP1-A_C1.tsv", sep = "\t", header=T)
msgfplus_output$index <- 1:nrow(msgfplus_output)

metadata <- read.table("data/mgf/mgf_metadata.txt") %>% as.matrix
begin <- which(metadata == "BEGIN")
end <- which(metadata == "END")

col_names <- c("TITLE", "RAWFILE", "SCANS", "CHARGE", "PEPMASS", "RTINSECONDS")
data <- matrix(ncol = length(col_names))[NULL,]
for (i in 1:length(begin)) {
  current <- metadata[(begin[i]+1):(end[i]-1),]
  current <- current %>% lapply(function(x) strsplit(x, split = "=") %>% unlist) %>% do.call(rbind, .) 
  current_data <- current[,2]
  names(current_data) <- current[,1]
  data <- rbind(data, current_data[col_names])
  print(i)
}

data <- data %>% as.data.frame
data_join <- data %>% select(SCANS, RTINSECONDS)
colnames(data_join) <- c("ScanNum", "RT")
data_join$ScanNum <- data_join$ScanNum %>% as.character %>% as.integer
data_join$RT
msgfplus_output_rt <- left_join(msgfplus_output, data_join, by = "ScanNum")

best_matches <- msgfplus_output_rt %>% select(ScanNum, MSGFScore, Peptide, Protein, RT)

# a <- strsplit(x = best_matches$Protein %>% as.character, split = "\\|")
# best_matches$Protein <- lapply(a, function(x) x[2]) %>% unlist

best_matches <- best_matches %>% arrange(ScanNum, desc(MSGFScore)) %>% 
  group_by(ScanNum) %>% 
  slice(1)

best_matches$
best_matches$MSGFScore %>% hist

best_matches_filter <- best_matches %>% filter(MSGFScore > 0)

best_matches_filter %>% group_by(Peptide) %>% summarise(count = n()) %>% arrange(desc(count))

msgfplus_output %>% filter(ScanNum == 720) %>% arrange(desc(MSGFScore)) %>% head(1) %>% .$MSGFScore

# best_matches <- list()
# i <- 1
# for (sn in scan_num) {
#   best_matches[[i]] <- filter(msgfplus_output, ScanNum == sn) %>% arrange(-MSGFScore) %>% head(1)
#   i <- i + 1
# }


mtcars %>% 
  select(gear, wt) %>% 
  arrange(gear, desc(wt)) %>% 
  group_by(gear) %>% 
  slice(seq(n()*.2))
