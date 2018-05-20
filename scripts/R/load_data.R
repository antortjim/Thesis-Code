source("settings.R")

data_list <- list()
i <- 1
for (f in files) {
  data_list[[i]] <- read_ms(f)
  i <- i + 1
}

experimental_design <- read.table(file.path(data_dir, "..", "experimental_design.tsv"), sep = "\t", header=T, stringsAsFactors = F)
experimental_design <- experimental_design %>% mutate(file = 1:nrow(.), rep = factor(rep))

data <- do.call(rbind, data_list) %>%
  full_join(., experimental_design, by = "file") %>%
  arrange_(condition_fields)

write.table(x = data, file = file.path(data_dir, "dataset.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = data$name %>% unique, file = file.path(output_dir, "efam_ids.txt"), quote = F, row.names = F, col.names=F)

