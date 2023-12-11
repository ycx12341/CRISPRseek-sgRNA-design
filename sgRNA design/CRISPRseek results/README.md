## CRISPRseek results ## 

This sub-folder contains the .R script for sgRNA design and the outputs it generates. 

**CRISPRseek attempt intron DMRs tbd.R**: this .R script demonstrates an example of calling CRISPRseek to generate sgRNAs for the target regions stored in the binary file **seq_intron_DMRs_tbd_valid_length.rds**.

**seq_intron_DMRs_tbd_valid_length.rds**: this binary file stores the DNA sequences of the target regions we would like to design sgRNAs for. 

**sgRNAs_introns_DMRs_xxx.rds**: these binary files are the outputs of sgRNAs design for the corresponding target regions. Note that CRISPRseek generates more detailed outputs for each target region automatically, which are stored in separate folders. Here in this sub-folder, we omitted these folders that store the more detailed outputs.  
