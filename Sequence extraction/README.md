## DNA sequence extraction ## 

This folder illustrates an example of obtaining the DNA sequences for the target regions. 

**Intron DMRs xxx to be designed.rds**: these files store the information of the target regions in the form of GRanges object. 

**Fasta export.R**: this .R script reads the .rds files that store the information of the target regions, extract the DNA sequences of them, and then export these sequences to a binary file. 

**seq_intron_DMRs_tbd_valid_length.rds**: this binary file contains the DNA sequences of the target regions we want to design sgRNAs for, which are stored as character strings. 
