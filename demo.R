library(tidyverse)
library(Seurat)
library(R.utils)

dir.create("/Users/kerencheng/Desktop/demo")
setwd("/Users/kerencheng/Desktop/demo")

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109037/suppl/GSE109037%5FAdultHumanSpermatogenesis%2Dreps1%2D3%5Ffilteredmatrix%2Emtx%2Egz", destfile = "matrix.mtx.gz")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109037/suppl/GSE109037%5FAdultHumanSpermatogenesis%2Dreps1%2D3%5Ffilteredgenes%2Etsv%2Egz", destfile = "genes.tsv.gz")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109037/suppl/GSE109037%5FAdultHumanSpermatogenesis%2Dreps1%2D3%5Ffiltered%5Fbarcodes%2Etsv%2Egz", destfile = "barcodes.tsv.gz")
# use mac or PC open the barcodes gz files will have unextend error due to a mistake when uploading the files. sorry!
gunzip("matrix.mtx.gz", remove=FALSE)
gunzip("genes.tsv.gz", remove=FALSE)
gunzip("barcodes.tsv.gz", remove=FALSE)

# Load the spermatogenesis dataset
spermatogenesis.data <- Read10X(data.dir = "/Users/kerencheng/Desktop/demo")
spermatogenesis.data[1:5, 1:5]
# Initialize the Seurat object with the raw (non-normalized data).
spermatogenesis <- CreateSeuratObject(counts = spermatogenesis.data, project = "spermatogenesis", min.cells = 3, min.features = 200)
spermatogenesis

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
spermatogenesis[["percent.mt"]] <- PercentageFeatureSet(spermatogenesis, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(spermatogenesis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(spermatogenesis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(spermatogenesis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

spermatogenesis <- NormalizeData(spermatogenesis, normalization.method = "LogNormalize", scale.factor = 10000)
spermatogenesis <- NormalizeData(spermatogenesis)
spermatogenesis <- FindVariableFeatures(spermatogenesis, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(spermatogenesis), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(spermatogenesis)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
#
all.genes <- rownames(spermatogenesis)
spermatogenesis <- ScaleData(spermatogenesis, features = all.genes)
spermatogenesis <- RunPCA(spermatogenesis, features = VariableFeatures(object = spermatogenesis))
# Examine and visualize PCA results a few different ways
print(spermatogenesis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(spermatogenesis, dims = 1:2, reduction = "pca")

DimPlot(spermatogenesis, reduction = "pca")
DimHeatmap(spermatogenesis, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(spermatogenesis, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
spermatogenesis <- JackStraw(spermatogenesis, num.replicate = 100)
spermatogenesis <- ScoreJackStraw(spermatogenesis, dims = 1:20)

JackStrawPlot(spermatogenesis, dims = 1:15)

ElbowPlot(spermatogenesis)

spermatogenesis <- FindNeighbors(spermatogenesis, dims = 1:16)
spermatogenesis <- FindClusters(spermatogenesis, 
                                reduction.type = "pca",
                                resolution = c(0.4, 0.8, 1.2),
                                dims.use = 1:16,
                                save.SNN = TRUE)
# Look at cluster IDs of the first 5 cells
head(Idents(spermatogenesis), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
spermatogenesis <- RunUMAP(spermatogenesis, dims = 1:10)
spermatogenesis <- RunTSNE(spermatogenesis, dims = 1:16)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(spermatogenesis, reduction = "umap")
DimPlot(spermatogenesis, reduction = "tsne") + theme_bw()






















