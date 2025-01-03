
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
# Load counts data
data <- read.csv("txi_counts_data.csv", row.names = 1)

# Load metadata
meta <- read.csv("metadata.csv", row.names = 1)


#plotting the counts data (counts vs no. of genes)
ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#analyzing the mean vs variance relationship
mean_counts <- apply(data[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

#obs1: we can clearly see from the above plot that mean!= variance and so we can't model using poisson distribution
#obs2: varaince is higher than mean for higher mean_counts even though the spread of data is very less at this point
#obs3: the variability (spread) of variance_counts is very high for low mean expression.
#inference: Therefore, we use Negative bionomial dist. and not poisson to study DGE for RNA seq data.

#we will not check for the presence of samples in the correct order in both datasets
all(colnames(data$counts) %in% rownames(meta)) #checks for presence
all(colnames(data$counts) == rownames(meta)) #checks for ordering

#similar check for txi files:
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

#STEP 1: Creating the DESeq2 object

dds2 <- DESeqDataSetFromMatrix(data, colData = meta, design = ~ sampletype)

#note: The design formula specifies the column(s) in the metadata table and how they should be used in the analysis.
      #out sampletype col. has 3 factor levels which tells DESeq2 to evaluate DGE w.r.t these 3 levels

View(counts(dds2))

#we now perform size factor estimation:
dds2 <- estimateSizeFactors(dds2)
sizeFactors(dds2)

#the list of size factors:
# Irrel_kd_1 Irrel_kd_2 Irrel_kd_3 Mov10_kd_2 Mov10_kd_3 
# 1.1150371  0.9606366  0.7493552  1.5634128  0.9359082 
# Mov10_oe_1 Mov10_oe_2 Mov10_oe_3 
# 1.2257749  1.1406863  0.6541689 

#getting the normalized count matrix:
normalized_counts <- counts(dds2, normalized=TRUE)
view(normalized_counts)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)
#quote=F: Prevents unwanted double quotes of the data.
#col.names=NA: Ensures the column names are correctly written in table and first cell at top left is left empty

#NOTE: DESeq2 DOESN'T USE THESE NORMALIZED COUNT DATA DIRECTLY FOR COMPUTING DGE. 
        #Instead DESeq2 uses the raw counts and models the normalization inside the Generalized Linear Model (GLM).
        #so, we can't use these normalzied counts to compute DGE with negative bionomial dist. directly!

#STEP 2: Quality control
#what to QC for? 
# 1. to check which samples are similar to one another, 
# 2. is the expt. design satisfied, 
# 3. sources of variation

#Sample-level QC with PCA/clustering can help us find the similarity btw samples easily and also outliers, if any

#STEP 2.1: rLOG transform (only for visualization purposes)
## NOTES:

#before performing PCA/clustering we need to log-transform the counts  and acheive variance stabilization. We use rlog transform approach.
#these PCA/clustering methods give best results if data is not heterskedastic and variance is not constant across all levels of indp. variable (value of one variable changes, then spread of other variable also changes)
#if pca is performed without variance stabilization, then PCA is results skew towards strongly exp. gene alone (since their variace is high)
#but just normal log2 transform doesn't work, as after that genes with low counts start skewing results due to inherent noise and show DGE wrongly
#so regularized log is used,where ridge penalty shrinks the log transform values to gene's average across all samples
#normal log2 transform 

rld <- rlog(dds2, blind=TRUE)
#blind=TRUE so that this transform is unbiased and doesn't consider sample info

#rld also returns a DESeq2 object since all paramters like size factors that was used for calc log trasnform are stored in this object

#PLOT PCA
plotPCA(rld, intgroup="sampletype") #PlotPCA is deseq2 function

#By default plotPCA() uses the top 500 most variable genes.

#just exploring other PCs:
rld_mat <- assay(rld) #this assay method converts the rld object into a dataframe. assay is a part of Summarizedexperiment package in DESeq2
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))

#Clustering analysis: Hierarchical Clustering
#this clustering using pheatmap() needs to use dataframe as input and not rld object. So we use rld_mat generated above

#calc sample-wise correlation
rld_cor <- cor(rld_mat)

#plotting the heatmap and using meta data for annotation:
pheatmap(rld_cor, annotation = meta)

#observation: both PCA and clustering show quality clustering within sample groups (seperation between sample groups) and prove that the data is of good quality

# STEP 3: Perform DGE with DESeq2
# we will feed the raw counts for this analysis. The NB model will be fitted for the raw counts and statistical analysis will be done.
dds3 <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
dds3 <- DESeq(dds3)

#in our dataset, we have three factor levels: control, mov10 overexp. and mov10 knockdown. 
# we are only interested in control vs other 2 factors
resultsNames(dds3)
contrast <- list(resultsNames(dds3)[1], resultsNames(dds3)[2])
results(dds3, contrast = contrast)

#speicfiying the contrast clearly:
#SYNTAX: contrast <- c("condition", "level_to_compare", "base_level")
#base level is out control

contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
res_tableOE <- results(dds3, contrast=contrast_oe, alpha = 0.05)
class(res_tableOE)

#visualize results:
res_tableOE %>% 
  data.frame() %>% 
  View()

# Filter genes by zero expression
res_tableOE[which(res_tableOE$baseMean == 0),] %>% 
  data.frame() %>% 
  View()


res_tableOE_unshrunken <- res_tableOE
# Apply fold change shrinkage
res_tableOE <- lfcShrink(dds3, coef="sampletype_MOV10_overexpression_vs_control", type="apeglm")

#Step 4: Filtering genes that do not meet the criteria before DGE analysis
#some genes will have  "NA" values for p-val/p-adj cols. these genes are:
  #a. genes with zero counts in all samples
  #b. Genes with an extreme count outliers (a count value in a sample for a particular gene is too high or too low)
  #c. Genes with a low mean normalized counts

# 4.a) Filter genes with zero counts in all samples:
res_tableOE[which(res_tableOE$baseMean == 0),] %>% 
  data.frame() %>% 
  View()
#so we filter out the genes by using baseMean==0 argument

# 4.b) filter Genes with an extreme count outliers:

#genes that are outliers tend to have "NA" values for p-val with a non-zero baseMean
res_tableOE[which(is.na(res_tableOE$padj) & is.na (res_tableOE$pvalue)
                  & res_tableOE$baseMean > 0), ] %>% 
  data.frame() %>% 
  View()

# 4.c) filter genes with low mean counts:
res_tableOE[which(!is.na(res_tableOE$pvalue) & 
                    is.na(res_tableOE$padj) & 
                    res_tableOE$baseMean > 0),] %>% 
  data.frame() %>% 
  View()

#MA Plot:
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
plotMA(res_tableOE, ylim=c(-2,2))

#STEP 5: Summary of results
summary(res_tableOE, alpha = 0.05)

#STEP 6: Extracting DEGs::
padj.cutoff <- 0.05
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#filtering the table to only keep significant genes:
sigOE <- res_tableOE_tb %>%
  dplyr::filter(padj < padj.cutoff)
View(sigOE)

# STEP 7:Visualization of DEGs:
mov10_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

#adding gene symbols to the normalized counts table

normalized_counts <- counts(dds3, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

grch38annot <- tx2gene %>% 
  dplyr::select(ensgene, symbol) %>% 
  dplyr::distinct()
normalized_counts <- merge(normalized_counts, grch38annot, by.x="gene", by.y="ensgene")

normalized_counts <- normalized_counts %>%
  as_tibble()

normalized_counts 

#Plotting DEGs:
#1. plot single gene counts across all samples
# Find the Ensembl ID of MOV10
grch38annot[grch38annot$symbol == "MOV10", "ensgene"] 
# [1] "ENSG00000155363"

plotCounts(dds3, gene="ENSG00000155363", intgroup="sampletype") 

#2. plot expression of single gene using ggplot2:
d <- plotCounts(dds3, gene="ENSG00000155363", intgroup="sampletype", returnData=TRUE)
d %>% View()

ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("MOV10") +
  theme(plot.title = element_text(hjust = 0.5))

#HEATMAP for few samples:
### Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_OEsig <- normalized_counts[,c(1:4,7:9)] %>% 
  dplyr::filter(gene %in% sigOE$gene)  

#plotting heatmap btw control and Mov10 over-expressed samples

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_OEsig[2:7], 
        
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

#Volcano plot:
#Create a logic vector based column to in results data that indicates whether a gene is DEG or not using padj:
res_tableOE_tb <- res_tableOE_tb %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)
#a threshold_OE col with TRUE/FALSE values is now created

## Volcano plot
ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

#STEP 8: Annotation:
library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
ah
# Query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
#I have downloaded and used pre-made annotations from annotations_hub
annotations_ahb <- read.csv("annotations_ahb.csv")
#Step 9: Over representation analysis
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)

#for performing these analysis, let's remove genes with p-adj==NA
res_tableOE_tb_noNAs <- filter(res_tableOE_tb, padj != "NA" )

## Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_tableOE_tb_noNAs, annotations_ahb, by=c("gene"="gene_id")) 
#with this merger we obtain info on type of gene (gene biotype), entrez_id as well ensembl_id in our data

#To perform over-representation analysis we need a list of background and list of significant genes

#the list of significant genes is derived using p_adj filter and bg genes is the list of all genes w/o applying p_adj filter 

#bg genes: 
allOE_genes <- as.character(res_ids$gene)

#significant genes:
sigOE <- dplyr::filter(res_ids, padj < 0.05)


sigOE_genes <- as.character(sigOE$gene)

#Now we perform GO enrichment analysis:
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "clusterProfiler_Mov10oe.csv")

#visualize results:
dotplot(ego, showCategory=50) #top 30 results

#STEP 10: GSEA:
## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]

#extract all fold changes:
foldchanges <- res_entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrezid
foldchanges <- sort(foldchanges, decreasing = TRUE)
view(foldchanges)

#GSEA:
set.seed(123456)

#GSEA using KEGG gene sets
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

gseaKEGG_results <- gseaKEGG@result
write.csv(gseaKEGG_results, "gseaOE_kegg_results.csv", quote=F)
View(gseaKEGG_results)

## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03008')
