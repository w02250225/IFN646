## Gene-level differential expression analysis using DESeq2
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
#instalation
for (p in c("tidyverse", "RColorBrewer", "pheatmap", "BiocManager", "ggplot2",
            "ggrepel"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p)
for (p in c("DESeq2", "DEGreport"))
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p)

#Load data
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1)
meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)
### Check classes of the data we just brought in
class(meta)
class(data)
View(meta)
View(data)

### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

# View the original counts matrix
View(counts(dds))

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

#PCA
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
### Plot PCA
plotPCA(rld, intgroup="sampletype")

# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))

#Hierarchical Clustering

### Extract the rlog matrix from the object
rld_mat <- assay(rld)
## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
rld_cor <- cor(rld_mat) ## cor() is a base R function
head(rld_cor) ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
pheatmap(rld_cor)

heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10,
         fontsize_row = 10, height=20)
#MOV10 DE analysis: exploring the dispersion estimates and assessing model fit

## Plot dispersion estimates
dds <- DESeq(dds)
plotDispEsts(dds)

#Differential expression analysis with DESeq2: model fitting and hypothesis testing
#Shrunken log2 foldchanges (LFC)
#design formula
design <- ~ sex + age + treatment

design <- ~ sex + age + treatment + sex:treatment

design <- ~ age + treat_sex

#MOV10 DE analysis

## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~ sampletype)
## Run analysis
dds <- DESeq(dds)












