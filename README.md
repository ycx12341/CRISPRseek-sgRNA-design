# CRISPRseek-sgRNA-design
 
This repository contains the workflow of short guide RNA (sgRNAs) design and selection, which iw performed using the R package CRISPRseek[^1]. In summary, the entire process is as follows. 

1. Extract the DNA sequences of the target regions we would like to design sgRNAs for.

2. Call CRISPRseek to design the sgRNAs for the target regions based on their DNA sequences, and select the ones that have the highest predicted efficacy.

See the README file within each folder for more details. 

[^1]: Zhu LJ, Holmes BR, Aronin N, Brodsky MH. CRISPRseek: a bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. PLoS One. 2014 Sep 23;9(9):e108424. doi: 10.1371/journal.pone.0108424. PMID: 25247697; PMCID: PMC4172692.
