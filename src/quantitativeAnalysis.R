###########################################################
#                                                         #
#            Differential peptide expression              #
#                     27/01/2017                          #
#                      M A Komor                          #
#                                                         #
###########################################################

# Rscript quantitativeAnalysis.R MaxQuant_output_peptides Quantitative_Analysis_group1_column_nr Quantitative_Analysis_group2_column_nr Extract_output_prefix Quantitative_analysis_imputation

### load functions

# install required packages
list.of.packages <- c("limma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
} 

# load
library(limma)

# Normalise peptide intensities to mean of samples' medians
median_normalisation <- function(a){
  return(mean(apply(a,2,median, na.rm=TRUE))/apply(a, 2, median, na.rm=TRUE))
}

normalize_df <- function(df){
  
  pep_log <- log10(df)
  pep_log[pep_log == -Inf] <- NA
  
  return(sweep(pep_log, 2, median_normalisation(pep_log), '*'))
}

# missing value imputation
imputation <- function(df){
  
  global_mean <- min(df, na.rm=TRUE)
  global_sd <- mean( apply( df, 1, sd, na.rm = TRUE), na.rm = TRUE)
  
  df[is.na(df)]<-rnorm(length(df[is.na(df)]), mean = global_mean, sd = global_sd)
  
  return(df)
}

# differential peptide expression
diffPepExpr <- function(df, gr1, gr2){
  design<-data.frame(c(rep(1,length(gr1+gr2))), c(rep(0,length(gr1)), rep(1,length(gr2))))
  
  fit <- lmFit(df, design)
  e = eBayes(fit)
  tab <- topTable(e,sort="none",number=dim(df)[1], coef=2)
  
  # add p-values to the table
  df$limma.logFC<-tab$logFC
  df$limma.p.value <- tab$P.Value
  df$limma.adj.p.value <- p.adjust(tab$P.Value, method = "BH")
  
  return(df)
  
}


### execute

main <- function(args){
  
  if (length(args) != 5){
    print("Please check the config file for:")
    print("[MaxQuant output] peptides = ")
    print("[Quantitative Analysis] group1_column_nr = ")
    print("[Quantitative Analysis] group2_column_nr = ")
    print("[Extract] output_prefix = ")
    print("[Quantitative_analysis] imputation = ")
    stop("There are missing arguments.")
  }

  
  peps <- read.delim(args[1])
  group1 <- as.numeric(unlist(strsplit(args[2], ",")))
  group2 <- as.numeric(unlist(strsplit(args[3], ",")))
  
  # get peptide intensities
  rawInt <- peps[, c(group1, group2)]
  rownames(rawInt) <- peps$Sequence
  
  # normalise
  p_log_norm <- normalize_df(rawInt)
  
  nr_samples <- as.data.frame(apply(p_log_norm, 1, function(x){
    length(x[!is.na(x)])
  }))
  colnames(nr_samples) <- c("NumberOfSamples")
  
  # imputation
  if (args[5] == "yes"){
    p_log_norm <- imputation(p_log_norm)
  }
  
  # get peptide table
  pep_sv <- read.delim(paste(args[4], "variantPeptides.txt", sep="/"))
  peptides <- subset(p_log_norm, row.names(p_log_norm) %in% pep_sv$Sequence)
  
  # differential peptide expression by limma
  peptides <- diffPepExpr(peptides, group1, group2)
  
  # merge tables as output
  peptideTable <- merge(pep_sv, peptides, by.x = "Sequence", by.y = "row.names")
  peptideTable <- merge(peptideTable, nr_samples, by.x = "Sequence", by.y ="row.names")

  write.table(peptideTable, paste(args[4], "variantPeptidesQA.txt", sep="/"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(peptideTable[peptideTable$Incl_Excl == "incl", ], 
              paste(args[4], "variantPeptidesInclusionQA.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(peptideTable[peptideTable$Incl_Excl == "excl", ], 
              paste(args[4], "variantPeptidesExclusionQA.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(peptideTable[peptideTable$NonCanonical == FALSE, ], 
              paste(args[4], "variantPeptidesCanonicalQA.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(peptideTable[peptideTable$NonCanonical == TRUE, ], 
              paste(args[4], "variantPeptidesNonCanonicalQA.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
}



main(commandArgs(TRUE))
