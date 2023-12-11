# Fasta export.R
# Author: Yunchen Xiao

# This .R file extracts the DNA sequences of the regions that we want to 
# design short guide RNAs (sgRNAs) for. 

# Set the workspace and load the necessary packages 
rm(list = ls())
library(BSgenome.Hsapiens.UCSC.hg19)
library(readr)
library(tictoc)
library(dplyr)
library(rtracklayer)
library(seqinr)
library(Biostrings)

# Read in the .rds file that stores the intron DMRs (in the format of GRanges) 
# we would like to design sgRNAs for
intron.DMRs.patwise.hyper.tbd <- read_rds("Intron DMRs Patwise Hyper to be designed.rds")
intron.DMRs.AL2.hyper.tbd <- read_rds("Intron DMRs AL2 Hyper to be designed.rds")
intron.DMRs.AL2.hypo.tbd <- read_rds("Intron DMRs AL2 Hypo to be designed.rds")
intron.DMRs.tbd <- c(intron.DMRs.patwise.hyper.tbd,
                     intron.DMRs.AL2.hyper.tbd,
                     intron.DMRs.AL2.hypo.tbd)

# Check the intron DMRs that have length < 23
ind.intron.DMRs.leq.23 <- which(lengths(intron.DMRs.tbd) < 23)

# Manually remove the intron DMRs that have lengths less than 23 base pairs. 
intron.DMRs.tbd.valid.length <- intron.DMRs.tbd

# Get the sequences of the intron DMRs contained in the file, then store them in 
# a list. 
intron.DMRs.seq.list <- list()
for (i in 1:length(intron.DMRs.tbd.valid.length)) {
  intron.DMRs.temp <- intron.DMRs.tbd.valid.length[i]
  start.temp <- start(intron.DMRs.temp)
  end.temp <- end(intron.DMRs.temp)
  seq.temp <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, 
                     names = as.character(intron.DMRs.temp@seqnames),
                     start = start.temp,
                     end = end.temp))
  intron.DMRs.seq.list <- append(intron.DMRs.seq.list, seq.temp)
}
names(intron.DMRs.seq.list) <- intron.DMRs.tbd.valid.length$introns_DMRs_ind

# Write the 6 sequences in one single .rds file
write_rds(intron.DMRs.seq.list, "seq_intron_DMRs_tbd_valid_length.rds")




