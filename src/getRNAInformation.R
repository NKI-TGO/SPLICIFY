###########################################
#                                         #
#  add isoform information from RNA level #
#                                         #
#               27/01/2017                #
#               M A Komor                 #
#                                         #
###########################################

# Rscript getRNAInformation.R RMATS_o matsToFasta_input matsToFasta_comparisonName Extract_output_prefix matsToFasta_ifMXE

### load functions

# read rMATS output
readrMATS <- function(RMATS_o, event, matsToFasta_input){
  filename <- paste(RMATS_o, "/MATS_output/", event, ".MATS.", matsToFasta_input, ".txt", sep="")
  if (file.exists(filename)){
    return(read.delim(filename))
  } else {
    stop(paste("File ", filename, " does not exist. Please check the config file.", sep=""))
  }
}

# fill in the peptide table
fillInTable<-function(df){
  
  # parse 
  thisiso <- unlist(strsplit(as.character(df[, 2]), "\\$"))
  
  coor <- strsplit(thisiso, ";")
  info <- sapply( strsplit( sapply(coor, "[[", 1), "%"), "[[", 1)
  df$ID <- sapply(strsplit(info, prefix), "[[", 2)
  df$Event <- sapply( strsplit( sapply(coor, "[[", 1), "%" ), "[[", 4)
  df$Gene <- sapply( strsplit( sapply(coor, "[[", 1), "%" ), "[[", 2)
  
  df$coordinates <- paste(sapply(coor, "[[", 2), sapply(coor, "[[", 3), 
                          sapply(coor, "[[", 4), sapply(coor, "[[", 5), sep = ";")
  
  
  # read RNA-seq results
  
  for (i in 1:length( df[, 1])){
    x <- eval( as.symbol( df$Event[i]))
    x <- x[x$ID == df$ID[i], ]
    
    df$FDR[i] <- x$FDR
    df$InclLevel1[i] <- mean(as.numeric(unlist(strsplit(as.character(x$IncLevel1), ","))))
    df$InclLevel2[i] <- mean(as.numeric(unlist(strsplit(as.character(x$IncLevel2), ","))))
    df$ExclLevel1[i] <- 1-df$InclLevel1[i]
    df$ExclLevel2[i] <- 1-df$InclLevel2[i]
    df$IncLevelDifference[i] <- x$IncLevelDifference
  }
  
  return(df)
  
}

### data analysis

# args <- c("/home/NFS/research_projects/tumor_specific_biomarkers/triplo/matssiSF3B1_rMATS3.2.5/", 
#          "ReadsOnTargetAndJunctionCounts", "siSF3B1", "data/prot/siSF3B1/", "yes")

# read arguments
main <- function(args){
  
  if (length(args) != 5){
    print("Please check the config file for:")
    print("[RMATS] o = ")
    print("[matsToFasta] input = ")
    print("[matsToFasta] comparisonName = ")
    print("[Extract] output_prefix = ")
    print("[matsToFasta] ifMXE = ")
    stop("There are missing arguments.")
  }
  
  prefix <<- args[3]
  
  # read rMATS splice variant information
  se <<- readrMATS(args[1], "SE", args[2])
  a3ss <<- readrMATS(args[1], "A3SS", args[2])
  a5ss <<- readrMATS(args[1], "A5SS", args[2])
  ri <<- readrMATS(args[1], "RI", args[2])
  if (args[5] == "yes" ) mxe <<- readrMATS(args[1], "MXE", args[2])
  
  # read PeptideTable
  peptideTable <- read.delim(paste(args[4], "variantPeptides.txt", sep="/"))
  
  if ("ID" %in% colnames(peptideTable)){
    stop("--------------Your variant peptide table already has RNA information.--------------")
  } 
  # add RNA information to the peptide table
  peptideTable <- fillInTable(peptideTable)
  
  # reshuffle 
  peptideTable <- peptideTable[, c("Sequence", "ID", "Event", "Gene", "Incl_Excl", 
                                   "coordinates", "FDR", "InclLevel1", 
                                   "ExclLevel1", "InclLevel2", "ExclLevel2", 
                                   "IncLevelDifference", "Sample", "NonCanonical", 
                                   "SecondVariant", "Peptide", "NumberOfIsoforms", 
                                   "Left", "Right", "Location")]
  
  
  # write output 
  write.table(peptideTable, paste(args[4], "variantPeptides.txt", sep="/"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(peptideTable[peptideTable$Incl_Excl == "incl", ], 
              paste(args[4], "variantPeptidesInclusion.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(peptideTable[peptideTable$Incl_Excl == "excl", ], 
              paste(args[4], "variantPeptidesExclusion.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(peptideTable[peptideTable$NonCanonical == FALSE, ], 
              paste(args[4], "variantPeptidesCanonical.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(peptideTable[peptideTable$NonCanonical == TRUE, ], 
              paste(args[4], "variantPeptidesNonCanonical.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)

}



# execute 
main(commandArgs(TRUE))







