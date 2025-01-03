library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)

#if (!requireNamespace("BiocManager", quietly = TRUE)) {
  #install.packages("BiocManager")
}
#BiocManager::install("DEGreport")


samples <- list.files(path = "./data", full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./data/", "") %>% 
  str_replace(".salmon", "")

tx2gene <- read.delim("tx2gene_grch38_ens94.txt")
tx2gene %>% View()


txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM")
txi$counts %>% View()
data <- txi$counts %>% 
  round() %>% 
  data.frame()


sampletype <- factor(c(rep("control",3), rep("MOV10_knockdown", 2), rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Save counts data to a CSV file
write.csv(data, "txi_counts_data.csv", row.names = TRUE)

# Save metadata to a CSV file
write.csv(meta, "metadata.csv", row.names = TRUE)
