##################################
#                                #
# Script for 3-frame translation #
#                                #
#           18/01/2017           #
#                                #
##################################

# Rscript threeFrame.R in.fasta out.fasta

# check if packages installed
list.of.packages <- c("seqinr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")


# load packages
library(seqinr)

# read arguments
args <- commandArgs(TRUE)

# read input
database <- read.fasta(file=args[1])

# 3 frame translation
frame1 <- lapply(database, translate, frame = 0, sens = "F")
frame2 <- lapply(database, translate, frame = 1, sens = "F")
frame3 <- lapply(database, translate, frame = 2, sens = "F")

# adjust fasta headers
attr(frame1, "name") <- paste( attr(frame1, "name") , "_1", sep="")
attr(frame2, "name") <- paste( attr(frame2, "name") , "_2", sep="")
attr(frame3, "name") <- paste( attr(frame3, "name") , "_3", sep="")

# save output
write.fasta(frame1, attr(frame1, "name"), args[2], open = "w", nbchar = 60, as.string = FALSE)
write.fasta(frame2, attr(frame2, "name"), args[2], open = "a", nbchar = 60, as.string = FALSE)
write.fasta(frame3, attr(frame3, "name"), args[2], open = "a", nbchar = 60, as.string = FALSE)


