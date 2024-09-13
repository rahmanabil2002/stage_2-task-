###########################################
#### Stage 2 task - HackBio Internship #### 
###########################################

# load the glioblastoma dataset
gbm <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv", row.names = 'X')

# explore dataset
head(gbm)
dim(gbm)
str(gbm)
table(is.na(gbm)) # no NAs

# 1.generate a heatmap for the entire dataset ####
library(gplots)
library(RColorBrewer) # for color palettes 

# diverging color palette -> for normalized data
div_color_palette <- rev(brewer.pal(11, "RdBu")) # red for up-regulation and blue for down-regulation

# sequential color palette -> for raw counts
seq_color_palette <- brewer.pal(9, "Blues") # blue gradient where light blue represents low counts

# distribution of raw counts across samples
boxplot(gbm, xlab = "samples", ylab = "counts", las = 2, col = "lightblue") # las = 2 is to rotate the x-axis labels 
hist(gbm[, "TCGA.19.4065.02A.11R.2005.01"], 
     main = "distribution of Raw Counts for Sample 1", 
     xlab = "counts", 
     col = "lightgreen", 
     breaks = 50)

# normalized data to reduce skewness of the data 
log_gbm <- log2(gbm + 1)
boxplot(log_gbm, xlab = "samples", ylab = "counts", las = 2, col = "lightblue") # we see the distribution clearer and eliminate the effect of outliners
hist(log_gbm[, "TCGA.19.4065.02A.11R.2005.01"], 
     main = "distribution of Raw Counts for Sample 1", 
     xlab = "log2(counts)", 
     col = "lightgreen", 
     breaks = 50)

# standardized data after log transformation to standardize each gene's expression across genes 
zscores_gbm <- t(scale(t(log_gbm))) # standardizes by rows (genes) so each gene gets a mean of 0 because I want to compare gene expression patterns across samples
boxplot(zscores_gbm, xlab = "samples", ylab = "counts", las = 2, col = "lightblue")
hist(zscores_gbm[, "TCGA.19.4065.02A.11R.2005.01"], 
     main = "distribution of Raw Counts for Sample 1", 
     xlab = "zscores", 
     col = "lightgreen", 
     breaks = 50)

# heatmap with sequential palette 
heatmap.2(as.matrix(gbm),
          col = seq_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          sepcolor = "black",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

heatmap.2(as.matrix(log_gbm),
          col = seq_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          sepcolor = "black",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

heatmap.2(as.matrix(zscores_gbm),
          col = seq_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          sepcolor = "black",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

# heatmap with diverging palette with standardized data
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

heatmap.2(as.matrix(log_gbm),
          col = div_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

heatmap.2(as.matrix(zscores_gbm),
          col = div_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

# heatmap clustering by genes
heatmap.2(as.matrix(zscores_gbm),
          col = div_color_palette,
          Rowv = T, Colv = F, dendrogram = "row",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

# heatmap clustering by samples
heatmap.2(as.matrix(zscores_gbm),
          col = div_color_palette,
          Rowv = F, Colv = T, dendrogram = "column",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

# heatmap clustering by genes and samples
heatmap.2(as.matrix(zscores_gbm),
          col = div_color_palette,
          Rowv = T, Colv = T, dendrogram = "both",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(10,10))

# 2.subset genes that are significantly upregulated and downregulated ####
# first 5 samples (02) are recurrent tumors, last 5 samples (01) are primary tumors
# PCA plot 
# change names of the samples in a duplicated df
gbm_copy <- zscores_gbm
colnames(gbm_copy) <- c("02A_1","02A_2","02A_3","02A_4","02A_5","01A_1","01A_2","01B_3","01A_4","01A_5")

pca <- prcomp(t(gbm_copy), center = TRUE, scale. = FALSE)
plot(pca) # how much variation is explained by each component 
summary(pca) 
pca$sdev # sd of PCs
PC1_and_PC2 <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], type = rownames(pca$x))
library(ggplot2)
ggplot(PC1_and_PC2,
       aes(x = PC1, y = PC2, col = type)) +
  geom_point() +
  geom_text(aes(label = type),
            hjust = 0, vjust = 0)

# create the df for "metadata"
sample_code <- c("TCGA.19.4065.02A.11R.2005.01", "TCGA.19.0957.02A.11R.2005.01", "TCGA.06.0152.02A.01R.2005.01", "TCGA.14.1402.02A.01R.2005.01", "TCGA.14.0736.02A.01R.2005.01", "TCGA.06.5410.01A.01R.1849.01", "TCGA.19.5960.01A.11R.1850.01", "TCGA.14.0781.01B.01R.1849.01", "TCGA.02.2483.01A.01R.1849.01", "TCGA.06.2570.01A.01R.1849.01")
tumor_stage <- c("recurrent", "recurrent", "recurrent", "recurrent", "recurrent", "primary", "primary", "primary", "primary", "primary")
meta <- data.frame(sample_code, tumor_stage)
meta
meta$tumor_stage <- as.factor(meta$tumor_stage)

# differential expression analysis
library(DESeq2)

dds <- DESeqDataSetFromMatrix( # dds is the DESeq2 dataset object
  countData = gbm,
  colData = meta,
  design = ~ tumor_stage
)

dds <- DESeq(dds) # differential expression analysis

res <- results(dds) # results of the differential expression analysis 
res

# filtering
up_genes <- res[which(res$padj < 0.05 & res$log2FoldChange > 1), ] # 2 times change in expression
up_genes@rownames # list of upregulated genes
down_genes <- res[which(res$padj < 0.05 & res$log2FoldChange < -1), ]
down_genes@rownames # list of downregulated genes 

background <- as.factor(res@rownames) # list of background genes

# Load necessary libraries
library(ggplot2)

# Load the data from the CSV file
data <- read.csv("D:/HackBio internship/Task stage 2/enrichment (1).csv")

# Display the first few rows to check the structure
head(data)
# Sort by Enrichment FDR and select top 5 pathways
top5_pathways <- data[order(data$Enrichment.FDR), ][1:5, ]

# Convert FDR to -log10(FDR)
top5_pathways$NegLogFDR <- -log10(top5_pathways$Enrichment.FDR)
# Create the lollipop plot
ggplot(top5_pathways, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
  geom_segment(aes(x = Pathway, xend = Pathway, y = 0, yend = nGenes), color = "grey") +
  geom_point(aes(size = NegLogFDR), color = "blue") +
  coord_flip() +
  labs(x = "Pathways", y = "Number of Genes", size = "-log10(FDR)") +
  theme_minimal()

