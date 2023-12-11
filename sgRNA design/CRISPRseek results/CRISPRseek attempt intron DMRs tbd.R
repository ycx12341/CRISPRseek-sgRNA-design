# CRISPRseek attempt intron DMRs tbd.R
# Author: Yunchen Xiao

# This .R file calls CRISPRseek and designs the sgRNAs for the targeting 
# sequences provided

# Set the workspace and load the necessary packages
rm(list = ls())
library(BSgenome.Hsapiens.UCSC.hg19)
library(readr)
library(CRISPRseek)
library(tictoc)
library(dplyr)
library(rtracklayer)
library(doParallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Load the .rds file that contains the 6 sequences to be investigated
intron.DMRs.seq.tbd <- read_rds("seq_intron_DMRs_tbd_valid_length.rds")

# Set up the cluster for parallel run. 
n.thread <- detectCores() - 1
n.sims <- length(intron.DMRs.seq.tbd)
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Set a directory to save the automatic results
# tic()
foreach(i = 1:n.sims, .packages = c("CRISPRseek", "BSgenome.Hsapiens.UCSC.hg19", 
                                    "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                                    "org.Hs.eg.db")) %dopar% {
  save.sims.dir <- paste0("Results_seek_m2_", names(intron.DMRs.seq.tbd)[i])
  if (dir.exists(save.sims.dir) == FALSE) {
    dir.create(save.sims.dir)
  }
  
  sgRNA.intron.DMRs.tbd.thres.10.m2.temp <- 
    offTargetAnalysis(inputFilePath = DNAStringSet(as.character(intron.DMRs.seq.tbd[i])),
                      BSgenomeName = BSgenome.Hsapiens.UCSC.hg19,
                      findgRNAs = TRUE,
                      findgRNAsWithREcutOnly = FALSE,
                      findPairedgRNAOnly = FALSE,
                      chromToSearch = paste0("chr", c(1:22, "X", "Y")),
                      gRNAoutputName = paste0("gRNA_seek_", names(intron.DMRs.seq.tbd)[i]),
                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                      outputDir = paste0("./", save.sims.dir),
                      max.mismatch = 2,
                      overwrite = TRUE,
                      orgAnn = org.Hs.egSYMBOL,
                      # n.cores.max = 1,
                      annotateExon = TRUE)
  
  readr::write_rds(sgRNA.intron.DMRs.tbd.thres.10.m2.temp, paste0("sgRNAs_", names(intron.DMRs.seq.tbd)[i], ".rds"))
}
toc()