read_msqrob <- function(filename) {
  RSqM_signif <- read.table(file = filename, header = T, sep = "\t", stringsAsFactors = F)
  return(RSqM_signif)
}
read_bayesquant <- function(filename) {
  # filename <- "../thp1/bayesquant_res.tsv"
  bayesquant <- read.table(filename, sep = "\t", header = T)
  colnames(bayesquant)[1] <- "protein"
  return(bayesquant)
}

read_quant <- function(filename, filetype="MSqRob") {
  if(filetype == "MSqRob") {
    quant <- read_msqrob(filename)
  } else if(filetype == "BayesQuant") {
    quant <- read_bayesquant(filename)
  }
  quant <- quant[!(1:nrow(quant) %in% grep(pattern = "CONTAMINANT", x = quant$Protein.IDs)),]
  return(quant)
}
