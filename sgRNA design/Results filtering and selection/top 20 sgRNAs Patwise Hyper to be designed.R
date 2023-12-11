# Intron DMRs Hyper Patwise to be designed top 20 sgRNAs
# Author: Yunchen Xiao

# Clear the current workspace and load the necessary packages
rm(list = ls())
library(readr)
library(GenomicRanges)
library(stringr)
library(dplyr)

# Read the intron DMRs (pat-spec hyper) to be designed
intron.DMRs.pat.hyper.to.be.designed <- 
  read_rds("Intron DMRs Patwise Hyper to be designed.rds")

# 1 of them!

# Check if the files containing the sequences we wanted to design sgRNAs for 
# are all included in the current directory
names.all.files <- list.files()
names.rds.files <- names.all.files[grep("sgRNAs_introns_DMRs_", names.all.files)]
names.rds.files.trun <- str_replace_all(names.rds.files, pattern = c("sgRNAs_|.rds"),
                                        replacement = "")
all.equal(sort(match(names.rds.files.trun, intron.DMRs.pat.hyper.to.be.designed$introns_DMRs_ind)),
          c(1:length(intron.DMRs.pat.hyper.to.be.designed))) # TRUE

# Filter the sgRNAs to keep the ones in the corresponding ranges only
for (i in 1:length(intron.DMRs.pat.hyper.to.be.designed)) {
  # Read in the corresponding intron DMRs
  intron.DMRs.temp <- intron.DMRs.pat.hyper.to.be.designed[i]
  
  # Read in the corresponding sgRNA information list
  sgRNAs.info.list.temp <- read_rds(paste0("sgRNAs_", intron.DMRs.temp$introns_DMRs_ind, ".rds"))
  
  # If valid sgRNAs were identified within the corresponding intron DMRs
  if (length(sgRNAs.info.list.temp) != 0) {
    # Extract the on.target sgRNAs
    on.target.temp <- sgRNAs.info.list.temp$on.target
    
    # Find the sgRNAs that are in the same chromosome as the intron DMR
    ind.same.chrom <- grep(paste0(as.character(seqnames(intron.DMRs.temp)), ":"), on.target.temp[,2])
    
    # If there are still valid sgRNAs left
    if (length(ind.same.chrom) != 0) {
      # Be careful when there is only one sgRNA left after the filter!
      if (length(ind.same.chrom) == 1) {
        on.target.same.chrom <- t(as.matrix(on.target.temp[ind.same.chrom,]))
      } else {
        on.target.same.chrom <- as.matrix(on.target.temp[ind.same.chrom,])
      }
      
      # Now extract those sgRNAs in the same chromosome and fall in the same 
      # range as the corresponding intron DMR
      on.target.same.chrom.valid <- vector()
      
      # For all the sgRNAs left in the same chromosome
      for (j in 1:length(on.target.same.chrom[,1])) {
        # Split the position information to see if it falls into the same range 
        # as the corresponding intron DMR
        on.target.same.chrom.temp <- on.target.same.chrom[j, ]
        range.split.temp <- unname(unlist(strsplit(on.target.same.chrom.temp[2], ":|-")))
        start.temp <- as.numeric(range.split.temp[2])
        end.temp <- as.numeric(range.split.temp[3])
        
        # If so, append the corresponding row into the matrix that stores the 
        # sgRNAs that are in the same chromosome and range as the corresponding 
        # intron DMR 
        if (start(intron.DMRs.temp) <= start.temp & end.temp <= end(intron.DMRs.temp)) {
          on.target.same.chrom.valid <- rbind(on.target.same.chrom.valid, on.target.same.chrom[j, ])
        }
      }
      
      # Append the matrix that stores the "valid" sgRNAs into the information list
      sgRNAs.info.list.temp[[(length(sgRNAs.info.list.temp) + 1)]] <- on.target.same.chrom.valid
      names(sgRNAs.info.list.temp)[length(sgRNAs.info.list.temp)] <- "on.target.valid"
      
      # Write the new .rds file
      write_rds(sgRNAs.info.list.temp, paste0("sgRNAs_", intron.DMRs.temp$introns_DMRs_ind, "_filtered.rds"))
    } else {
      # If there are no sgRNAs in the same chromosome as the corresponding intron DMR, 
      # set the filtered information list to be empty and write it down as an empty .rds file
      sgRNAs.info.list.temp.valid <- vector()
      write_rds(sgRNAs.info.list.temp.valid, paste0("sgRNAs_", intron.DMRs.temp$introns_DMRs_ind, "_filtered.rds"))
    } 
  } else {
    # If there are no sgRNAs in identified in the corresponding intron DMR, 
    # set the filtered information list to be empty and write it down as an empty .rds file
    sgRNAs.info.list.temp.valid <- vector()
    write_rds(sgRNAs.info.list.temp.valid, paste0("sgRNAs_", intron.DMRs.temp$introns_DMRs_ind, "_filtered.rds"))
  }
  # Optional, can be used to track the progress
  print(i)
}

# Now get the top 20 sgRNAs for these filtered results
top.20.sgRNAs.list <- vector(mode = "list", length = length(intron.DMRs.pat.hyper.to.be.designed))

for (i in 1:length(top.20.sgRNAs.list)) {
  # Read in the current .rds file that stores the sgRNAs identified for the intron DMR
  res.temp <- read_rds(paste0("sgRNAs_", intron.DMRs.pat.hyper.to.be.designed[i]$introns_DMRs_ind, "_filtered.rds"))
  # If there are valid sgRNAs identified
  if (length(res.temp) != 0) {
    # Extract the on and off-target information of the sgRNAs
    sgRNAs.on.tar.temp <- res.temp$on.target.valid
    if(dim(sgRNAs.on.tar.temp)[1] > 1) {
      sgRNAs.on.tar.temp <- unique(sgRNAs.on.tar.temp[,2:4])
    }
    sgRNAs.off.tar.temp <- res.temp$offtarget
    # If there are more than 20 valid sgRNAs identified
    if (length(sgRNAs.on.tar.temp[, 1]) > 20) {
      # Sort the top 20 sgRNAs
      ind.top.20.on.tar.temp <- order(sgRNAs.on.tar.temp[,3], decreasing = TRUE)[1:20]
      
      # Extract the top 20 sgRNAs from the on target matrix
      sgRNAs.on.tar.top.20.temp <- sgRNAs.on.tar.temp[ind.top.20.on.tar.temp,]
      
      # Extract the full information of these top 20 sgRNAs from the off target matrix
      sgRNAs.on.tar.top.20.full.temp <- left_join(as.data.frame(sgRNAs.on.tar.top.20.temp),
                                                  as.data.frame(sgRNAs.off.tar.temp),
                                                  by = c("forViewInUCSC" = "forViewInUCSC",
                                                         "extendedSequence" = "extendedSequence"))
      
      # Generate a GRanges object to store these sgRNAs
      sgRNAs.on.tar.top.20.GRanges.temp <- GRanges(seqnames = sgRNAs.on.tar.top.20.full.temp$chrom,
                                                   ranges = IRanges(start = as.numeric(sgRNAs.on.tar.top.20.full.temp$chromStart),
                                                                    end = as.numeric(sgRNAs.on.tar.top.20.full.temp$chromEnd)),
                                                   names = sgRNAs.on.tar.top.20.full.temp$name,
                                                   strand = sgRNAs.on.tar.top.20.full.temp$strand,
                                                   score = sgRNAs.on.tar.top.20.full.temp$score,
                                                   efficacy = sgRNAs.on.tar.top.20.full.temp$gRNAefficacy.y,
                                                   n.mismatch = sgRNAs.on.tar.top.20.full.temp$n.mismatch,
                                                   gRNA.seq = sgRNAs.on.tar.top.20.full.temp$gRNAPlusPAM)
      top.20.sgRNAs.list[i] <- unique(sgRNAs.on.tar.top.20.GRanges.temp)
      print(i)
    } else {
      # If there are less than 20 sgRNAs identified in the current intron DMR, 
      # then extract all their information
      sgRNAs.on.tar.full.temp <- left_join(as.data.frame(sgRNAs.on.tar.temp),
                                           as.data.frame(sgRNAs.off.tar.temp),
                                           by = c("forViewInUCSC" = "forViewInUCSC",
                                                  "extendedSequence" = "extendedSequence"))  
      
      # Store these information in a GRange list
      sgRNAs.on.tar.GRanges.temp <- GRanges(seqnames = sgRNAs.on.tar.full.temp$chrom,
                                            ranges = IRanges(start = as.numeric(sgRNAs.on.tar.full.temp$chromStart),
                                                             end = as.numeric(sgRNAs.on.tar.full.temp$chromEnd)),
                                            names = sgRNAs.on.tar.full.temp$name,
                                            strand = sgRNAs.on.tar.full.temp$strand,
                                            score = sgRNAs.on.tar.full.temp$score,
                                            efficacy = sgRNAs.on.tar.full.temp$gRNAefficacy.y,
                                            n.mismatch = sgRNAs.on.tar.full.temp$n.mismatch,
                                            gRNA.seq = sgRNAs.on.tar.full.temp$gRNAPlusPAM)
      top.20.sgRNAs.list[i] <- unique(sgRNAs.on.tar.GRanges.temp)
      print(i)
    }
  } else {
    print("No valid sgRNAs identified in the current intron DMR")
    top.20.sgRNAs.list[i] <- "No valid sgRNAs identified in the current intron DMR"
    print(i)
  }
}
names(top.20.sgRNAs.list) <- paste0("top_20_sgRNAs_", intron.DMRs.pat.hyper.to.be.designed$introns_DMRs_ind)

# Export the top 20 sgRNAs selected into a binary file
write_rds(top.20.sgRNAs.list, "Top 20 sgRNAs intron DMRs patwise hyper incomplete.rds")
