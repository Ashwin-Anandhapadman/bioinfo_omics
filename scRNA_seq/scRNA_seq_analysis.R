# Single-cell RNA-seq analysis - QC
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(SingleCellExperiment)

#Step0: just a basic brush-up
ctrl_counts <- Read10X(data.dir = "data") #Loading just control data
ctrl_counts[1:3,1:5]
print(dim(ctrl_counts))
#33538 genes 737280 cells

ctrl <- CreateSeuratObject(counts = ctrl_counts, min.features = 100)
# View metadata
head(ctrl@meta.data)
slotNames(ctrl@assays$RNA)

#checking data features
print(dim(ctrl)) #33538 genes 15688 cells
slotNames(ctrl)
#misc: Miscellaneous additional data.
#layers: Likely contains data layers such as counts, normalized counts, etc.

#metadata:
#1. nCount_RNA: Number of unique molecular identifiers (UMIs) per cell.
#2. nFeature_RNA: Number of genes detected per cell.

#trying to load only .tsv files gives error:
ctrl_ash <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")
ctrl_ash_seu <- CreateSeuratObject(counts = ctrl_counts, min.features = 100)
# View metadata
head(ctrl_ash_seu@meta.data)
#the above code gives error: 
# Error in Read10X(data.dir = "my_data_c") : Barcode file missing. Expecting barcodes.tsv.gz



#step1: creating the Seurat object for both sample files 

for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}
#we have min.features as 100 becuz we want to filter out cells with less than 100 genes expressed

print(dim(seurat_obj)) #33538 15756

#Step2: checking the metadata in the new Seurat object
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)

print(dim(ctrl_raw_feature_bc_matrix)) #33538 15688
print(dim(stim_raw_feature_bc_matrix)) #33538 15756

#Step3: merging both seurat objects (sample and control)
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, 
                       y = stim_raw_feature_bc_matrix, 
                       add.cell.id = c("ctrl", "stim"))
# Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat) 

head(merged_seurat@meta.data) #ctrl data
tail(merged_seurat@meta.data) #sample data

print(dim(merged_seurat)) #33538 31444 (row stacking of ctrl and sample)

#step 4: QC
View(merged_seurat@meta.data)

#finding the no. of genes per UMI and adding the info to Seurat meta data
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
#info on how calc is done in PercentageFeatureSet:
    #For each cell, the function takes the sum of counts across all mt- genes (features) and then   
    #divides by the count sum for all genes (incld mt genes for each cell). This value is multiplied by 100 to obtain a percentage value.We don't need percentage, so we divide by 100 to get ratio

merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

#step5:Create a seperate metadata dataframe
metadata <- merged_seurat@meta.data

# a new cols of cell name (called "cells") is added. it is obtained by using the row indices of the dataframe itself
metadata$cells <- rownames(metadata)

#creating a new column called sample that defines whether the sample is ctrl or patient sample
metadata$sample <- NA #create an empty col
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add and storing the metadata back to Seurat object
merged_seurat@meta.data <- metadata
save(merged_seurat, file="data/merged_filtered_seurat.RData")

#Step 5: Assessing quality metrics:
#1. Cell counts
#2. UMI counts per cell
#3. Genes detected per cell
#4. Complexity (novelty score)
#5. Mitochondrial counts ratio

# 5.1. Cell counts: Ideally we expect # of unique barcodes to be equal to #cells in sample
#not always the case because more than one barcode can be inside a bead or there can be barcodes without cells or dying cells too

#Visualize no. of cell counts per sample (x-axis is the new sample cols. we made):
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

#observation1: we see slightly more than 15k cells in both samples and this is more than the expected cell count (~13k).
#observation2: therefore, there could be presence of some junk cells in the samples that need to be handled.


#5.2. UMI counts per cell:

ggplot(metadata, aes(color=sample, x=nUMI, fill= sample)) + 
geom_density(alpha = 0.3) + 
scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
geom_vline(xintercept = 500)
#alpha is just a measure of transparency in plot 

#note: we need to ideally have more than 500 UMI counts per cell
#obs1: we can see from the density plot that most of our cells have more than 1000 UMIs, which is a good sign.

#5.3. Genes per cell:
  ggplot(metadata, aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
  
#obs: we have good quality cell data as most cells have more than 300 genes per cell
  
#Genes per sample:
  metadata %>% 
    ggplot(aes(x=sample,y=nGene, color=sample, fill=sample)) + 
    geom_density(alpha=0.2) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NGene")
  
#5.4. Novelty score:
  #note: novelty score is the ratio btw nGene to UMI counts indicating complexity. 
          #cells with very low gene counts but high UMI counts are low complex cells becuz they just happen to have high UMI due to seq. depth
          # we have already calc. this ratio and stored in our metadata as log10GenesPerUMI
  
metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)

#obs: we can see our cells have high complexity based on the expected novelty score thereshold of 0.80
  
#5.5 Mito counts ratio:


metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#obs: Based on the expected mito ratio threshold of 0.2 we can see our cells are of high quality with very low mito ratio estimates

#Note: Always consider joint filtering effects. QC metrics shd not be considered in isolation always!

  ggplot(metadata, aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
  
  #obs: this plot clearly compared the gene counts against UMI counts against the backdrop of mito ratio.
  #obs2: we can clearly see a good positive correlation between UMI count and gene count especially in the upper right quadrant, indicating high data quality
  #obs3: the lower left quadrant mostly represents the low quality cells with relatively high mito contamination too.
  
#step 6: Filtering the samples based on the metrics

  #6.1 cell-level filtering

  #based on the widely followed filtering metrics, we can use the below thresholds:
  
  #nUMI > 500
  #nGene > 250
  #log10GenesPerUMI > 0.8
  #mitoRatio < 0.2
  
#we use the subset function to achieve this filtering
filtered_seurat <- subset(x = merged_seurat, 
                            subset= (nUMI >= 500) & 
                              (nGene >= 250) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.20))
  
  
  #6.2. Gene-level filtering:
#to remove genes with zero counts'

# Extract the counts separately
counts_gene <- GetAssayData(object = filtered_seurat, layer = "counts")
nonzero <- counts > 0

#keep only genes which are expressed in 10 or more cells.
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

#final filtered seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

save(filtered_seurat, file="seurat_filtered_new.RData")
save(filtered_seurat, file="~/seurat_filtered.RData")

#Re-analysis for checking:
metadata_clean <- filtered_seurat@meta.data

#UMI counts per cell:

ggplot(metadata_clean, aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.3) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

#mito ratio:

metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#step 7: Normalization
seurat_phase <- NormalizeData(filtered_seurat)
load("data/cycle.rda")

# we want to analyze the effect of cell cycle on transcription as it is a source of variation
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
View(seurat_phase@meta.data)

#To assign each cell a score based on its expression of G2/M and S phase markers,
#we can use the Seurat function CellCycleScoring()

#After scoring for cell cycle, we use PCA to find out whethe cell cycle is a key driver of variation or not
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     verbose = FALSE)

seurat_phase <- ScaleData(seurat_phase)

#identifying 15 most highly variable genes:
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]
print(top_genes)
#list of top 15 genes:

# "HBB"    "HBA2"   "CCL4L2" "CCL7"   "IGKC"   "HBA1"  
# "PPBP"   "CCL4"   "CCL3"   "CCL8"   "TXN"    "CXCL10"
# "CXCL1"  "CCL2"   "IFI27" 



#plotting the variable genes
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE, xnudge=0, ynudge=0)


seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")

#We do not see large differences due to cell cycle phase. 
#Based on this plot, we would not regress out the variation due to cell cycle.

#Step 8: Clustering


DimHeatmap(seurat_phase, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = seurat_phase[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = seurat_phase, 
          ndims = 40)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_phase, 
                                   dims = 1:40)

seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

seurat_integrated@meta.data %>% 
View()

Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Compute UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

# Plot UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        group.by = "seurat_clusters")

table(Idents(seurat_integrated))
head(seurat_integrated@meta.data)


markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, 
                          min.pct = 0.1, logfc.threshold = 0.25)
