library("dplyr")
library("readxl")
library("dplyr")
library("tidyr")
library("ggplot2")
theme_set(theme_bw())
library("viridis")
library("tibble")

base_dir <- ifelse(
  Sys.info()[["sysname"]] == "Linux",
  "/z/home/aoj/thesis/genedata",
  "//hest/aoj/thesis/genedata/")

# project_name <- "Saphera"
project_name <- "thp1"
# compare <- "commercialvspurified"
# analysis_conditions <- list(condition = c("commercial", "purified"))
analysis_conditions <- list(condition1 = c("supernatant"))
# condition_fields <- c("condition")
condition_fields <- c("condition1", "condition2")
# database_accession <- "generic"
database_accession <- "uniprot"
  
project_dir <- file.path(base_dir, project_name)
data_dir <- file.path(project_dir, "data", "xlsx")
plots_dir <- file.path(project_dir, "plots")
output_dir <- file.path(project_dir, "output")

ggsave2 <- function(filename, plot) {
  print(filename)
  return(ggsave(filename = file.path(plots_dir, filename), plot = plot))
}
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

read_ms <- function(filename) {
  # Reads the excel file in filename
  file_path <- file.path(data_dir, filename)
  print(filename)
  dat <- read_excel(file_path)
  metadata_fields <- c("Source Type", "Group", "Description", "Maximum RT Shift", "Flags", "Search Results", "Enzyme", "Fraction")
  metadata_rows <- which(dat$Name %in% metadata_fields)
  metadata <- dat[metadata_rows, 1:2]
  dat <- dat[-metadata_rows,]
  
  dat <- dat[-(dat$Name %>% grep(pattern = "con_| : ", x = .)), ]
  
  colnames(dat)[2] <- "Amount"
  colnames(dat) <- colnames(dat) %>% make.names() %>% tolower
  numeric_cols <- apply(dat, 2, function(x) as.numeric(x) %>% is.na %>% all %>% !.) %>% which
  dat[,numeric_cols] <- apply(dat[,numeric_cols], 2, as.numeric)
  dat$file <- i
  dat <- dat[,c("name", "amount", "file", "protein.p.value")]
  return(dat)
}



excluded_proteins <- "P33BJB"
dep_criteria <- c("BOTH")
palette <- rep("red", 4)


log.fold.change.threshold <- 1
p.value.threshold <- -log10(0.05)
correlation.threshold <- 0.96

files <- file.path(data_dir) %>% list.files %>% grep(pattern = "^\\w.*\\.xlsx", value = T)

