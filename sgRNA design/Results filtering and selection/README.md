## Results filtering and selection ##

This sub-folder contains an example of filtering and selecting the top sgRNAs generated for the corresponding target regions. 

**Intron DMRs Patwise Hyper to be designed.rds**: this binary file contains the GRanges object for the target region(s). It is presented in the sub-folder to ensure sgRNAs are designed for all related target regions.

**top 20 sgRNAs Patwise Hyper to be designed**: this .R script demonstrates an example of filtering the sgRNAs designed by CRISPRseek, keep the ones within the corresponding target regions and select the top ones based on the predicted efficacy. 

**sgRNAs_intron_DMRs_xxx.rds**: these binary files store the sgRNAs designed by CRISPRseek initially and the valid ones within the corresponding target regions after filtering. 

**Top 20 sgRNAs intron DMRs patwise hyper incomplete.rds**: this binary file stores the top 20 sgRNAs selected from the filtered valids ones (**sgRNAs_intron_DMRs_patwise_hyper_46_filtered.rds**). 

