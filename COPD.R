library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratWrappers)
library(scater)
library(Matrix)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(Rtsne)
library(cowplot)
library(ggplot2)
library(ggsci)
library(scales)
library(MAST)
library(DOSE)
library(patchwork)
library(plotly)
library(monocle)
library(MASS)
library(loomR)
library(RColorBrewer)
library(grDevices)
library(colorRamps)
library(data.table)
library(hexbin)


### Import 10x dataset as a seurat object


data_dir <- '~/Analysis/COPD/JK01'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix1 <- Read10X(data.dir = data_dir)
jk01 = CreateSeuratObject(counts = expression_matrix1, project = "JK01")

data_dir <- '~/Analysis/COPD/JK02'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix2 <- Read10X(data.dir = data_dir)
jk02 = CreateSeuratObject(counts = expression_matrix2, project = "JK02")

data_dir <- '~/Analysis/COPD/JK03'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix3 <- Read10X(data.dir = data_dir)
jk03 = CreateSeuratObject(counts = expression_matrix3, project = "JK03")

data_dir <- '~/Analysis/COPD/JK04'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix4 <- Read10X(data.dir = data_dir)
jk04 = CreateSeuratObject(counts = expression_matrix4, project = "JK04")

data_dir <- '~/Analysis/COPD/JK05'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix5 <- Read10X(data.dir = data_dir)
jk05 = CreateSeuratObject(counts = expression_matrix5, project = "JK05")

data_dir <- '~/Analysis/COPD/JK06'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix6 <- Read10X(data.dir = data_dir)
jk06 = CreateSeuratObject(counts = expression_matrix6, project = "JK06")

data_dir <- '~/Analysis/COPD/JK07'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix7 <- Read10X(data.dir = data_dir)
jk07 = CreateSeuratObject(counts = expression_matrix7, project = "JK07")

data_dir <- '~/Analysis/COPD/JK08'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix8 <- Read10X(data.dir = data_dir)
jk08 = CreateSeuratObject(counts = expression_matrix8, project = "JK08")

data_dir <- '~/Analysis/COPD/JK09'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix9 <- Read10X(data.dir = data_dir)
jk09 = CreateSeuratObject(counts = expression_matrix9, project = "JK09")

data_dir <- '~/Analysis/COPD/JK10'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix10 <- Read10X(data.dir = data_dir)
jk10 = CreateSeuratObject(counts = expression_matrix10, project = "JK10")

data_dir <- '~/Analysis/COPD/JK11'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix11 <- Read10X(data.dir = data_dir)
jk11 = CreateSeuratObject(counts = expression_matrix11, project = "JK11")

data_dir <- '~/Analysis/COPD/JK12'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix12 <- Read10X(data.dir = data_dir)
jk12 = CreateSeuratObject(counts = expression_matrix12, project = "JK12")



COPDmeta <- read_csv("COPDmeta.csv")
COPDmeta = data.frame(COPDmeta)

md = NULL
meta = NULL
metaall = list()
JK = c("jk01", "jk02", "jk03", "jk04", "jk05", "jk06", "jk07", "jk08", "jk09", "jk10", "jk11", "jk12")

for(i in 1:length(JK)){
	
		md = unlist(COPDmeta[i,])
		s <- paste(sprintf(JK[i]))

		meta = eval(parse(text = s))@meta.data
		meta$Class <- md[3]
		meta$Age <- md[4]
		meta$Gender <- md[5]
		meta$BMI <- md[6]
		meta$BrinkmanIndex <- md[7]
		meta$FEV1.FVC <- md[8]
		meta$FEV1 <- md[9]
		meta$FEV1.pred <- md[10]
		meta$Cormobidities <- md[11]

	metaall[i] = list(meta)
	}

jk01 = CreateSeuratObject(counts = expression_matrix1, project = "JK01", meta.data = metaall[[1]])
jk02 = CreateSeuratObject(counts = expression_matrix2, project = "JK02", meta.data = metaall[[2]])
jk03 = CreateSeuratObject(counts = expression_matrix3, project = "JK03", meta.data = metaall[[3]])
jk04 = CreateSeuratObject(counts = expression_matrix4, project = "JK04", meta.data = metaall[[4]])
jk05 = CreateSeuratObject(counts = expression_matrix5, project = "JK05", meta.data = metaall[[5]])
jk06 = CreateSeuratObject(counts = expression_matrix6, project = "JK06", meta.data = metaall[[6]])
jk07 = CreateSeuratObject(counts = expression_matrix7, project = "JK07", meta.data = metaall[[7]])
jk08 = CreateSeuratObject(counts = expression_matrix8, project = "JK08", meta.data = metaall[[8]])
jk09 = CreateSeuratObject(counts = expression_matrix9, project = "JK09", meta.data = metaall[[9]])
jk10 = CreateSeuratObject(counts = expression_matrix10, project = "JK10", meta.data = metaall[[10]])
jk11 = CreateSeuratObject(counts = expression_matrix11, project = "JK11", meta.data = metaall[[11]])
jk12 = CreateSeuratObject(counts = expression_matrix12, project = "JK12", meta.data = metaall[[12]])



COPD <- merge(jk01, y = c(jk02, jk03, jk04, jk05, jk06, jk07, jk08, jk09, jk10, jk11, jk12), project = "COPD")


### Filtering 
COPD[["percent.mt"]] <- PercentageFeatureSet(COPD, pattern = "^MT-")

VlnPlot(COPD, features = c("nFeature_RNA"))
VlnPlot(COPD, features = c("nCount_RNA"))
VlnPlot(COPD, features = c("percent.mt"))


plot1 <- FeatureScatter(COPD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(COPD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

COPD <- subset(COPD, subset = nFeature_RNA > 100 & percent.mt < 75)

COPD <- NormalizeData(COPD, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(COPD)
COPD <- ScaleData(COPD, vars.to.regress = "percent.mt")

COPD <- FindVariableFeatures(COPD)

top10 <- head(VariableFeatures(COPD), 10)
plot1 <- VariableFeaturePlot(COPD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2




s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

COPD <- CellCycleScoring(COPD, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



### PC plot
COPD <- RunPCA(COPD, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

DimPlot(COPD, reduction = "pca")
DimPlot(COPD, reduction = "pca", group.by = "Class")

DimHeatmap(COPD, dims = 1:15, cells = 500, balanced = TRUE)

VizDimLoadings(COPD, dims = 1:2, reduction = "pca")


### Jackstraw plot
COPD <- JackStraw(COPD, num.replicate = 100)
COPD <- ScoreJackStraw(COPD, dims = 1:20)

JackStrawPlot(COPD, dims = 1:15)

ElbowPlot(COPD)



### find Cluster
COPD <- FindNeighbors(COPD, dims = 1:15)
COPD <- FindClusters(COPD, resolution = 0.2)

### tSNE
COPD <- RunTSNE(COPD, dims = 1:75)

DimPlot(COPD, label = TRUE, label.size = 10, reduction = "tsne") + NoLegend()

colpall = primary.colors(17)
colpall = sample(colpall)
DimPlot(COPD, label = TRUE, label.size = 8, cols = colpall, reduction = "tsne") + NoLegend()


colclass = c("red3", "deepskyblue3", "chocolate1")

tp1 = DimPlot(COPD, reduction = "tsne") 
tp2 = DimPlot(COPD, reduction = "tsne", group.by = "Class", cols = colclass)
plot_grid(tp1, tp2)

tp3 = DimPlot(COPD, reduction = "tsne", label = TRUE, label.size = 8) + NoLegend()
tp4 = DimPlot(COPD, reduction = "tsne", group.by = "Phase")
plot_grid(tp3, tp4)

DimPlot(COPD, reduction = "tsne", group.by = "Class", cols = colclass)


### save dataset
saveRDS(COPD, file = "COPD.rds")




#### Density plot
COPD.Dim = DimPlot(COPD, reduction = "tsne")
mat = data.frame(UMAP1 = COPD.Dim$data$tSNE_1, UMAP2 = COPD.Dim$data$tSNE_2, ident = COPD.Dim$data$ident)

bin　<-　hexbin(COPD.Dim$data$tSNE_1, COPD.Dim$data$tSNE_2, xbins = 100)
mCols = colorRampPalette(rev(brewer.pal(64, 'Spectral')))
plot(bin, main = "", colramp = mCols, legend = FALSE) 


rd = kde2d(mat[,1], mat[,2], n = 2700)
image(rd)
contour(rd)




