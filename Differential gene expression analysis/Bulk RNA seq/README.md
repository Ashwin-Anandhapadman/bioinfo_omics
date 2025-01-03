# Differential gene expression analysis with bulk RNA seq data using DESeq2

This is my project on implementation of differential gene expression (DGE) analysis on bulk RNA seq data. This project is a part of my industrial training in bioinformatics.
The dataset used in this project is based on this paper: http://www.ncbi.nlm.nih.gov/pubmed/25464849. 

In this paper, the RNA-Seq analysis was performed on HEK293F cells subjected to three distinct experimental conditions:
1. Mov10 oe: Cells transfected with a MOV10 transgene (overexpression).
2. Mov10 kd: Cells treated with siRNA targeting Mov10 to knock down its expression.
3. Irrelevant kd: Cells treated with non-specific siRNA as a control.

The authors of the paper investigated the interactions among various genes associated with Fragile X syndrome (FXS), a genetic disorder characterized by the abnormal production of the fragile X mental retardation protein (FMRP). MOV10 protein is a RNA helicase that is associated with FMRP.The dataset includes multiple replicates for each condition. 

This project aims to explore the transcriptional patterns associated with the modulation of MOV10 expression, with the irrelevant siRNA condition serving as our control.

## The following steps were performed in this project:
1. DESeq2 object creation and QC analysis on count data using PCA and clustering
2. Differential gene expression analysis using DESeq2
3. Visualization of differentially expressed genes.
4. Functional analysis (Over-representation analysis and Gene set enrichment analysis (GSEA))

