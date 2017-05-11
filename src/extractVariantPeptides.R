#################################################
#                                               #
# Extract variant peptides from MaxQuant output #
#                 25/01/2017                    #
#                 M A Komor                     #
#                                               #
#################################################

# Rscript extractVariantPeptides.R evidence.txt database.pep.fasta path/to/output/ processor_nr canonicalDB.fasta 

# install required packages
list.of.packages <- c("foreach", "doParallel", "seqinr", "cwhmisc", "reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load libraries
library(foreach)
library(doParallel)
library(seqinr)
library(cwhmisc)
library(reshape)

### load functions

# fasta database list to data.frame
getFastaDf <- function(databaseList){
  databaseDf <- foreach(i = 1:length(databaseList), .combine = rbind, .packages=c('seqinr') ) %dopar% {
    name <- attr( databaseList[i], 'name')
    seq <- toupper( paste( getSequence( databaseList[i])[[1]], collapse = ''))
    data.frame( name, seq)
  }
  return(databaseDf)
}

# get peptides that are in sequenced from the splice variants database
getPeptidesFromVariants <- function(allpeptides, databaseDf){
  foundPeptides <- foreach (i = 1:length(allpeptides[, 1]), .combine = rbind) %dopar% {
    
    ident <- databaseDf[ grep( allpeptides[i, 1], databaseDf[, 2]), 1]
    seq <- allpeptides[i, 1]
    gene <- paste( ident, collapse = '$')
    
    # return peptide and fasta headers
    if ( length(ident) > 0 ){
      data.frame( seq, gene)
    }  
  }
  return(foundPeptides)
  
}

# Evaluate if peptide is unique for the splice variant
isUnique <- function(found){
  foundUnique <- foreach (i = 1:length(found[,1]), .combine = rbind, .packages = c('reshape')) %dopar% {
    
    # uniq <- 0 as default means not unique (shared) peptide
    uniq <- 0
    
    # get all the splice variants the peptide maps to
    ident1 <- unlist( strsplit( as.character( found[i,2]), '\\$'))
    ident <- sapply( strsplit( ident1, ';'), '[[', 1)
    ids_event <- paste( sapply(strsplit(ident, '%'), '[[', 1), 
                        sapply(strsplit(ident, '%'), '[[', 4), sep = '%')
    
    seq <- found[i, 1]
    gene <- found[i, 2]
    
    # form splice variant data.frame for the peptide
    all <- data.frame(id = NA, gene = NA, which = NA, event = NA )
    for (j in 1:length( ident )){
      all[j,] <- unlist( strsplit( ident[j], "\\%"))
    }
    all$id <- paste( all$id, all$event, sep = '%')
    
    # data.frame for each splive variant event in a row, gives incl and excl value
    # indicating if peptide occurs in incl, excl variant. If incl = 0, 
    # peptide is specific for excl, if excl = 0, peptide is specific for incl
    casted <- cast(all, id~which, fun.aggregate = length, value = 'event')
    
    if ( dim( casted)[2] == 3){
      uniq_iso <- na.omit( data.frame( id = NA, incl = NA, excl = NA))
      
      # check if peptide specific for any of the isoforms 
      if ( length(casted[casted$excl == 0,]) != 0) uniq_iso <- do.call('rbind', list(uniq_iso, casted[casted$excl == 0,]))
      if ( length(casted[casted$incl == 0,]) != 0) uniq_iso <- do.call('rbind', list(uniq_iso, casted[casted$incl == 0,]))
      
      if (dim(uniq_iso)[1] != 0){
        
        # uniq <- 1 means peptide is unique
        uniq <- 1
        
        # get fasta header for which peptide is specific
        for (id_nr in 1:length(uniq_iso$id)){
          if (id_nr == 1){
            gene <- ident1[ grep( uniq_iso$id[ id_nr], ids_event)]
          }else{
            gene <- paste ( gene, ident1[ grep (uniq_iso$id[id_nr], ids_event)], sep = '$')
          }
        }
      }  
      
      # if there are no incl or excl columns (so only gene and incl or excl column)
      # it means that the peptide is specific for the incl/excl form 
    }else if (dim(casted)[2] < 3) uniq <- 1
    
    # return data.frame with uniq column incidating if peptide is splice variant specific
    data.frame(seq, gene, uniq)
  }
  
  return(foundUnique)
}

# get sample information from evidence.txt file
getSamples <- function(peptideDf, evidence){
  outputDf <- foreach (i=1:length(peptideDf[,1]), .combine = rbind) %dopar% {
    
    # in which Experiment was the peptide found
    sample<-paste( unique( evidence[ 
      grep( peptideDf[i,1], evidence[,1]), 'Experiment'])[ order(
        unique (evidence[ grep( peptideDf[i,1], evidence[,1]),'Experiment']))], 
      collapse=';')
    
    # return data.frame
    data.frame(peptideDf[i,], sample)
  }
  return(outputDf)
}

# get uniq_pep table with separate rows for each isoform
# adds additional rows if a peptide maps to multiple isoforms

dfTransform <- function(uniq_pep){
  newdf<-NULL
  
  for (i in 1:dim( uniq_pep)[1]){
    div <- unique( unlist( strsplit( as.character( uniq_pep[i,2]), '\\$')))
    
    for (k in div){
      newdf <- rbind( newdf, data.frame( as.character( uniq_pep[i,1]), k, uniq_pep[i, 3:dim( uniq_pep)[2]], row.names = NULL))
    }
  }
  
  colnames(newdf)<-colnames(uniq_pep)
  return(newdf)
}

# get location of each peptide
peptideLocation <- function(newdf, database){
  for (row in 1:dim(newdf)[1]){
  
    is<-newdf[row,]
    
    # get splice variant sequence
    db_row <- database[as.character(database[,1])==as.character(is[,2]),]
    
    loc <- unlist ( strsplit ( as.character(is$gene) , ';' ) )[ 4 ]
    index <- c ( ) 
    
    # get frame of translation
    frame <- as.numeric ( unlist ( strsplit ( unlist ( strsplit ( as.character(is$gene) , ';' ) )[ 5  ] , '_' )) [ 2 ] )
    
    for ( i in unlist ( strsplit ( loc, '@') ) ) {
      
      # segment start and end coordinates
      startend <- unlist ( strsplit ( i , '-' ) )
      index <- append ( index,  seq ( ( 1 + as.numeric ( startend[ 1 ] ) ), as.numeric ( startend[ 2 ] ) , 1) )
    }
    
    if (unlist ( strsplit ( as.character(is$gene) , ';' ) )[ 3 ] == '-') index <- rev ( index )
    
    pep_ind <- index [ seq ( frame, length ( index ) , 3 ) ]
    
    start <- cpos ( as.character(db_row[, 'seq']) , as.character(is[ , 'seq' ]), 1)
    
    inds <- pep_ind [ start:( start -1 + nchar ( as.character ( is [ , 'seq' ] ) ) ) ]
    
    
    left<-1
    right<-1
    if (unlist ( strsplit ( as.character(is$gene) , ';' ) )[ 3 ] =='+'){
      for (i in  1:(length(inds)-1)){
        if ((inds[i+1]-inds[i])==3) left <- left + 1
        else break
      }
      
      for (i in  ((length(inds))):2){
        if ((inds[i]-inds[i-1])==3) right <- right + 1
        else break
      }
    }else{
      for (i in  1:(length(inds)-1)){
        if ((inds[i+1]-inds[i])==-3) left <- left + 1
        else break
      }
      
      for (i in  ((length(inds))):2){
        if ((inds[i]-inds[i-1])==-3) right <- right + 1
        else break
      }
    }
    
    newdf[row, 'left']<-left
    newdf[row, 'right']<-right
    newdf[row, 'loc']<-paste(inds, collapse=';')
    
  }
  return(newdf)
}

# is peptide split, spanning or on target
peptideAnnotate <- function(newdf){
  for (i in 1:dim(newdf)[1]){
    
    
    if ( nchar( as.character( newdf[i,'seq'])) >= 
         ( as.numeric( newdf[i,'left']) + as.numeric( newdf[i,'right']))){
      
      newdf[i, 'pep']<-'split peptide'
      
    }else if ( grepl( "incl", newdf[i, 2] )){
      
      newdf[i, 'pep']<-'on target'
      
      loc <- unlist ( strsplit ( as.character( newdf[i,2] ), ';' ))[ 4 ]
      mid <- unlist ( strsplit ( as.character( loc), "@" )) [2]
      start <- as.numeric( unlist ( strsplit (mid, "-"))[1] )
      end <- as.numeric( unlist ( strsplit (mid, "-"))[2] )
      
      split_loc <- as.character( c( start-1, start, start+1, end-1, end, end+1 ))
      inside_loc <- unlist( strsplit( newdf[i, "loc"], ";")) [2:( 
        length( unlist( strsplit( newdf[i, "loc"], ";" ))) - 1)] 
      
      for (j in split_loc){
        if (j %in% inside_loc) newdf[i, 'pep']<-'spanning peptide'
        
      }
    }else{
      newdf[i, 'pep']<-'on target'
    }
  }
  return(newdf)
}

# find second variant, excl for incl, and incl for excl
getSecondVariant <- function(x, df){
  
  isoform<-sapply(strsplit(as.character(x[2]), ";"), "[[", 1)
  iso.vec <- unlist(strsplit(isoform, "%"))
  
  opposite <- ifelse(iso.vec[3] == "incl", "excl", "incl")
  second.var <- paste(iso.vec[1], iso.vec[2], opposite, iso.vec[4], sep="%")
  
  return(TRUE %in% grepl(second.var, df[,2]))
  
}


##### main

main <- function(args){
  
  if (length(args) != 5){
    print("Please check the config file for:")
    print("[MaxQuant output] evidence = ")
    print("[matsToFasta] outputPepFasta = ")
    print("[Extract] output_prefix = ")
    print("[Extract] threads = ")
    print("[Extract] canonical = ")
    stop("There are missing arguments.")
    
  }
  
  ### read parameters
  # read canonical sequence database
  canonical <- read.fasta(file=args[5])
  
  # load database with splice variants
  database <- read.fasta(file=args[2])
  
  # load MaxQuant evidence.txt file
  ev <- read.delim(args[1], header=TRUE)
  
  # If experiment not defined in MaxQuant, MS raw file name used as experiment
  if (length (grep('Experiment', colnames(ev))) == 0 ) ev$Experiment <- ev$Raw.file
  
  ### setup parallel backend to use X processors
  cl <- makeCluster(as.numeric(args[4]))
  registerDoParallel(cl)
  
  ### data analysis
  
  # process the fasta file
  database <- getFastaDf(database)
  canonical <- getFastaDf(canonical)
  
  # get unique peptide sequences
  peps <- data.frame(seq=unique(ev[,1]))
  
  # get peptides from splice variant database
  svPeptides <- getPeptidesFromVariants(peps, database)
  
  # remove peptide list to save on memory
  rm(peps)
  
  # Evaluate if peptide is unique for the splice variant
  svPeptides <- isUnique(svPeptides)
  
  # get only specific peptides
  uniq_pep <- svPeptides[svPeptides[,3] == 1,]
  rm(svPeptides)
  
  # get sample information from evidence.txt file
  uniq_pep <- getSamples(uniq_pep, ev)
  rm(ev)
  
  # identify reference peptides
  rp <- getPeptidesFromVariants(uniq_pep, canonical)
  rm(canonical)
  
  # stop parallel computing
  stopCluster(cl)
  
  # tranform the dataframe
  peptideTable <- dfTransform(uniq_pep)
  
  # add peptide location
  peptideTable <- peptideLocation(peptideTable, database)
  
  # peptide information
  peptideTable <- peptideAnnotate(peptideTable)
  
  # remove column uniq
  peptideTable <- peptideTable[, -3]
  
  # identify NonCanonical peptides
  noncanpep <- subset(peptideTable, !(peptideTable[,1] %in% rp[,1]))
  canpep <- subset(peptideTable, (peptideTable[,1] %in% rp[,1]))
  
  peptideTable$NonCanonical <- ifelse(peptideTable[,1] %in% noncanpep[,1], TRUE, FALSE)
  
  # identify incl and excl peptides
  incl <- peptideTable[grepl("incl", peptideTable[,2]),]
  excl <- peptideTable[grepl("excl", peptideTable[,2]),]
  
  peptideTable$Incl_Excl <- ifelse(grepl("incl", peptideTable[,2]), "incl", "excl")
  
  peptideTable$NumbrOfIsoforms <- apply(peptideTable, 1, function(x){
    return(length(peptideTable[peptideTable[,1] == x[1],1]))
  })
  
  peptideTable[,2] <- as.character(peptideTable[,2])                           
  peptideTable$SecondVariant <- apply(peptideTable, 1, getSecondVariant, df = peptideTable) 
  
  colnames(peptideTable) <- c("Sequence", "SpliceVariant", "Sample", "Left", 
                              "Right", "Location", "Peptide", "NonCanonical", "Incl_Excl", 
                              "NumberOfIsoforms", "SecondVariant")
  
  
  
  write.table(peptideTable, paste(args[3], "variantPeptides.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(peptideTable[peptideTable$Incl_Excl == "incl", ], 
              paste(args[3], "variantPeptidesInclusion.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(peptideTable[peptideTable$Incl_Excl == "excl", ], 
              paste(args[3], "variantPeptidesExclusion.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(peptideTable[peptideTable$NonCanonical == FALSE, ], 
              paste(args[3], "variantPeptidesCanonical.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(peptideTable[peptideTable$NonCanonical == TRUE, ], 
              paste(args[3], "variantPeptidesNonCanonical.txt", sep = "/"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
}

### execute
main(commandArgs(TRUE))



