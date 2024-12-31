# scRNA analysis 

This my project on bioinformatic analysis of single cell RNA seq data.

The dataset is derived from a large study performed by Kang et al. (https://www.nature.com/articles/nbt.4042) whereby the authors present a computational algorithm that harnesses genetic variation (eQTL) to determine the genetic identity of each droplet containing a single cell (singlet) and identify droplets containing two cells from different individuals (doublets).The researchers tested their algorithm using pooled Peripheral Blood Mononuclear Cells (PBMCs) from eight lupus patients. These samples were divided into two groups: one control group and one group treated with interferon beta. The main aim of this study is to perform scRNA seq analysis on the stimulated (IFN-β) and control lupus patient samples and study patterns and variation in gene expression. The interferon treatment was performed in the original study to analyze the effect of the IFN-β treatment especially on immune cells in Lupus patients.

![image](https://github.com/user-attachments/assets/8bf6a1fa-f920-4bc2-b0a2-8c501801acab)

More info about the dataset: 

1. Libraries were created using the version 2 chemistry from 10X Genomics.
2. Sequencing was carried out on the Illumina NextSeq 500.
3. PBMC samples were collected from eight individual lupus patients and split into two separate aliquots.
4. One aliquot was treated with 100 U/mL of recombinant IFN-β for a duration of 6 hours.The second aliquot remained untreated.
5. After excluding doublets, 12,138 cells were identified in the control sample and 12,167 in the stimulated sample.\

 

The unprocessed dataset and the cell cycle data along was obtained from the Harvard Bioinformatics core's (HBC) tutorial on scRNA seq analysis: https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html. I would like to thank this HBC group for creating such a comprehensive and useful tutorial on performing scRNA analysis.
