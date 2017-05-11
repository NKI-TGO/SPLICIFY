##################################################
#                                                #
#                  MATS TO BED                   #
#  script to obtain bed file from rMATS output   #
#                                                #
#                  24/01/2017                    #    
#                   MA Komor                     #
#                                                #
##################################################

### How to run
# Rscript matsToBed.R path/to/rMATS/output comparisonName path/to/output.bed ifMXE input FDR
# e.g. Rscript matsToBed.R ../data/rmats/MATSoutput/ cancervsnormal ../data/rmats/output.bed yes ReadsOnTargetAndJunctionCounts 0.05

# Read arguments
args <- commandArgs(TRUE)


# Assign to variables
exp <- args[2]
fdr <- as.numeric(args[6])

# ### Skipped exon
se <- read.delim( paste (args[1], '/SE.MATS.', args[5], '.txt', sep = '') )
s <- se[se$FDR <= fdr,]

# data frame to bed formatting
if (nrow(s) != 0) {
	# inclusion variants
	bed <- data.frame(chr=s$chr)

	bed[,2] <- s$upstreamES
	bed[,3] <- s$downstreamEE
	bed[,4] <- paste( paste( paste(exp, s$ID, sep = ''),
	                         s$geneSymbol, 'incl', 'se', sep = '%'),
	                  s$chr, s$strand,
	                  paste( paste (s[,8],s[,9], sep = '-'),
	                         paste (s[,6],s[,7], sep = '-') ,
	                         paste (s[,10],s[,11], sep = '-'),
	                         sep = '@' ), '', sep = ';')
	bed[,5] <- 1
	bed[,6] <- s$strand
	bed[,7] <- s$upstreamES
	bed[,8] <- s$downstreamEE
	bed[,9] <- 0
	bed[,10] <- 3
	bed[,11] <- paste( s[,9] - s[,8], s[,7] - s[,6], s[,11] - s[,10], sep = ',')
	bed[,12] <- paste( 0,
	                   s$exonStart_0base - s$upstreamES,
	                   s$downstreamES - s$upstreamES,
	                   sep = ',')

	incl <- bed

	# exclusion variant
	bed <- data.frame(chr = s$chr)

	bed[,2] <- s$upstreamES
	bed[,3] <- s$downstreamEE
	bed[,4] <- paste ( paste(paste(exp, s$ID, sep = ''),
	                         s$geneSymbol, 'excl', 'se', sep = '%'),
	                   s$chr, s$strand,
	                   paste ( paste (s[,8],s[,9], sep = '-'),
	                           paste (s[,10],s[,11], sep = '-'),
	                           sep = '@'), '', sep = ';')
	bed[,5] <- 0
	bed[,6] <- s$strand
	bed[,7] <- s$upstreamES
	bed[,8] <- s$downstreamEE
	bed[,9] <- 0
	bed[,10] <- 2
	bed[,11] <- paste( s[,9] - s[,8], s[,11] - s[,10], sep = ',')
	bed[,12] <- paste( 0, s$downstreamES - s$upstreamES, sep = ',')


	sebed <- do.call(rbind, list(incl, bed))
}

### Retained intron
se <- read.delim( paste (args[1], '/RI.MATS.', args[5], '.txt', sep = '') )
s <- se[se$FDR<=fdr,]

# data frame to bed formatting
if (nrow(s) != 0){
  # inclusion variant
	bed <- data.frame(chr=s$chr)

	bed[,2] <- s$upstreamES
	bed[,3] <- s$downstreamEE
	bed[,4] <- paste ( paste ( paste( exp, s$ID, sep = ''),
	                           s$geneSymbol, 'incl', 'ri', sep = '%'),
	                   s$chr, s$strand ,
	                   paste(paste(s[,8], s[,9], sep = '-'),
	                         paste( (s[,9]) , s[,10] , sep = '-' ),
	                         paste( s[,10], s[,11], sep = '-') ,
	                         sep = "@" ),'', sep = ';')
	bed[,5] <- 1
	bed[,6] <- s$strand
	bed[,7] <- s$upstreamES
	bed[,8] <- s$downstreamEE
	bed[,9] <- 0
	bed[,10] <- 1
	bed[,11] <- s[,7] - s[,6]
	bed[,12] <- 0

	incl <- bed

	# exclusion variant
	bed <- data.frame(chr = s$chr)

	bed[,2] <- s$upstreamES
	bed[,3] <- s$downstreamEE
	bed[,4] <- paste ( paste(paste(exp, s$ID, sep = ''),
	                         s$geneSymbol, 'excl', 'ri', sep = '%'),
	                   s$chr, s$strand,
	                   paste ( paste (s[,8],s[,9], sep = '-'),
	                           paste (s[,10],s[,11], sep = '-') ,
	                           sep = '@' ), '', sep = ';')
	bed[,5] <- 0
	bed[,6] <- s$strand
	bed[,7] <- s$upstreamES
	bed[,8] <- s$downstreamEE
	bed[,9] <- 0
	bed[,10] <- 2
	bed[,11] <- paste(s[,9] - s[,8], s[,11] - s[,10], sep = ',')
	bed[,12] <- paste(0, s$downstreamES - s$upstreamES, sep = ',')

	sebed <- do.call(rbind, list(sebed,incl, bed))
}


se <- read.delim( paste (args[1], '/A5SS.MATS.', args[5], '.txt', sep = '') )
s <- se[se$FDR<=fdr,]

if (nrow(s) != 0){
  
  # inclusion variant
	bed <- data.frame(chr = s$chr)

	bed[,2] <- 1
	bed[,3] <- 1

	fasta_forward <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'incl', 'a5ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,8],s[,9], sep = '-'), 
	                                 paste (s[,9],s[,7], sep = '-'), 
	                                 paste (s[,10],s[,11], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	fasta_reverse <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'incl', 'a5ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,10],s[,11], sep = '-'), 
	                                 paste (s[,6],s[,8], sep = '-'), 
	                                 paste (s[, 8],s[, 9], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	bed[,4] <- ""
	bed[,5] <- 1
	bed[,6] <- s$strand

	
	bed[,4] <- ifelse(bed[,6] == "-", fasta_reverse, fasta_forward )
	
	
	bed[,2] <- ifelse(bed[,6] == '+', s$longExonStart_0base, s$flankingES)
	bed[,3] <- ifelse(bed[,6] == '+', s$flankingEE, s$longExonEnd)

	bed[,7] <- s$flankingES
	bed[,8] <- s$longExonEnd
	bed[,9] <- 0
	bed[,10] <- 2

	bed[,11] <- ifelse(bed[,6] == '+', 
	                   paste(s[,7] - s[,6], s[,11] - s[,10], sep = ','), 
	                   paste(s[,11] - s[,10], s[,7] - s[,6], sep = ','))
	bed[,12] <- paste(0, abs(s$flankingES - s$longExonStart_0base), sep = ',')

	incl<-bed

	# exclusion variant
	bed <- data.frame(chr = s$chr)
	bed[,2] <- 1
	bed[,3] <- 1

	fasta_forward <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'excl', 'a5ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,8],s[,9], sep = '-'), 
	                                 paste (s[,10],s[,11], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	fasta_reverse <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'excl', 'a5ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,10],s[,11], sep = '-'), 
	                                 paste (s[, 8],s[, 9], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	
	bed[,4] <-""
	
	bed[,5] <- 0
	bed[,6] <- s$strand
	
	bed[,4] <- ifelse(bed[,6] == "-", fasta_reverse, fasta_forward )
	
	bed[,2] <- ifelse(bed[,6] == '+', s$shortES, s$flankingES)
	bed[,3] <- ifelse(bed[,6] == '+', s$flankingEE, s$shortEE)

	bed[,7] <- min(s$shortES, s$flankingES)
	bed[,8] <- max(s$flankingEE, s$shortEE)
	bed[,9] <- 0
	bed[,10] <- 2

	bed[,11] <- ifelse(bed[,6] == '+', 
	                   paste(s[,9] - s[,8], s[,11] - s[,10], sep = ','), 
	                   paste(s[,11] - s[,10], s[,9] - s[,8], sep = ','))
	bed[,12] <- paste(0, abs(s$flankingES - s$shortES), sep = ',')

	sebed <- do.call(rbind, list(sebed,incl, bed))
}

se <- read.delim( paste (args[1], '/A3SS.MATS.', args[5], '.txt', sep = '') )
s <- se[se$FDR<=fdr,]

if (nrow(s) != 0){
  
  # inclusion variant
	bed <- data.frame(chr=s$chr)

	bed[,2] <- 1
	bed[,3] <- 1
	bed[,4] <- ""
	
	fasta_reverse <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'incl', 'a3ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,8],s[,9], sep = '-'), 
	                                 paste (s[,9],s[,7], sep = '-'), 
	                                 paste (s[,10],s[,11], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	fasta_forward <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'incl', 'a3ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,10],s[,11], sep = '-'), 
	                                 paste (s[,6],s[,8], sep = '-'), 
	                                 paste (s[, 8],s[, 9], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	
	bed[,5] <- 1
	bed[,6] <- s$strand
	
  bed[,4] <- ifelse(bed[,6] == "-", fasta_reverse, fasta_forward)
	
	bed[,2] <- ifelse(bed[,6] == '-', s$longExonStart_0base, s$flankingES)
	bed[,3] <- ifelse(bed[,6] == '-', s$flankingEE, s$longExonEnd)

	bed[,7] <- s$flankingES
	bed[,8] <- s$longExonEnd
	bed[,9] <- 0
	bed[,10] <- 2
	
	bed[,11] <- ifelse(bed[,6] == '-', 
	                   paste(s[,7] - s[,6], s[,11] - s[,10], sep = ','), 
	                   paste(s[,11]-s[,10], s[,7]-s[,6], sep = ','))
	bed[,12] <- paste(0, abs(s$flankingES - s$longExonStart_0base), sep = ',')

	incl <- bed

	# exclusion variant
	bed <- data.frame(chr = s$chr)
	bed[,2] <- 1
	bed[,3] <- 1

	fasta_reverse <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'excl', 'a3ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,8],s[,9], sep = '-'), 
	                                 paste (s[,10],s[,11], sep = '-'), 
	                                 sep = '@'), '', sep = ';')
	
	fasta_forward <- paste ( paste(paste(exp, s$ID, sep = ''), 
	                               s$geneSymbol, 'excl', 'a3ss', sep = '%'), 
	                         s$chr, s$strand, 
	                         paste ( paste (s[,10],s[,11], sep = '-'), 
	                                 paste (s[, 8],s[, 9], sep = '-'), 
	                                 sep = '@'), '', sep = ';')

	bed[,4] <- ""
	
	bed[,5] <- 0
	bed[,6] <- s$strand
	
	bed[,4] <- ifelse(bed[,6] == "-", fasta_reverse, fasta_forward)
	

	bed[,2] <- ifelse(bed[,6] == '-', s$shortES, s$flankingES)
	bed[,3] <- ifelse(bed[,6] == '-', s$flankingEE, s$shortEE)

	bed[,7] <- min(s$shortES, s$flankingES)
	bed[,8] <- max(s$flankingEE, s$shortEE)
	bed[,9] <- 0
	bed[,10] <- 2

	bed[,11] <- ifelse(bed[,6] == '-', 
	                   paste(s[,9] - s[,8], s[,11] - s[,10], sep = ','), 
	                   paste(s[,11] - s[,10], s[,9] - s[,8], sep = ','))
	bed[,12] <- paste(0, abs(s$flankingES - s$shortES), sep = ',')

	sebed <- do.call(rbind, list(sebed,incl, bed))
}

# Including mutually exclusive exons in the output
if (args[4] == 'yes'){

  se <- read.delim( paste (args[1], '/MXE.MATS.', args[5], '.txt', sep = '') )
  s <- se[se$FDR <= fdr,]

  if (nrow(s) != 0){
    #inclusion of exon 1 tagged as incl
    bed <- data.frame(chr = s$chr)

    bed[,2] <- s$upstreamES
    bed[,3] <- s$downstreamEE
    bed[,4] <- paste ( paste(paste(exp, s$ID, sep = ''),
                             s$geneSymbol, 'incl', 'mxe', sep = '%'),
                       s$chr, s$strand,
                       paste ( paste (s[,10],s[,11], sep = '-'),
                               paste (s[,6],s[,7], sep = '-'),
                               paste (s[,12],s[,13], sep = '-'),
                               sep = '@'), '', sep = ';')
    bed[,5] <- 1
    bed[,6] <- s$strand
    bed[,7] <- s$upstreamES
    bed[,8] <- s$downstreamEE
    bed[,9] <- 0
    bed[,10] <- 3
    bed[,11] <- paste(s[,11] - s[,10],
                      s[,7] - s[,6],
                      s[,13] - s[,12], sep = ',')
    bed[,12] <- paste(0,
                      s$X1stExonStart_0base - s$upstreamES,
                      s$downstreamES - s$upstreamES, sep = ',')

    incl <- bed

    #inclusion of exon 2 tagged as excl
    bed <- data.frame(chr = s$chr)

    bed[,2] <- s$upstreamES
    bed[,3] <- s$downstreamEE
    bed[,4] <- paste ( paste( paste(exp, s$ID, sep = ''),
                              s$geneSymbol, 'excl', 'mxe', sep = '%'),
                       s$chr, s$strand,
                       paste ( paste (s[,10],s[,11], sep = '-'),
                               paste (s[,8],s[,9], sep = '-'),
                               paste (s[,12],s[,13], sep = '-'),
                               sep = '@'), '', sep = ';')
    bed[,5] <- 0
    bed[,6] <- s$strand
    bed[,7] <- s$upstreamES
    bed[,8] <- s$downstreamEE
    bed[,9] <- 0
    bed[,10] <- 3
    bed[,11] <- paste(s[,11] - s[,10],
                      s[,9] - s[,8],
                      s[,13] - s[,12], sep = ',')
    bed[,12] <- paste(0, s$X2ndExonStart_0base - s$upstreamES,
                      s$downstreamES - s$upstreamES, sep = ',')

    sebed <- do.call(rbind, list(sebed,incl, bed))
  }
}

# save the results
write.table(sebed , args[3] , row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

