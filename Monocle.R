### Monocle

### Ready to perform monocle
### Extraction of Alveolar subpopulations
sel = c(0, 1, 2, 3, 5, 8, 9, 14, 15, 19)
ALV = subset(EPI, ident = sel)
saveRDS(ALV, file = "ALV.rds")


### Extracton of Bronchiolar subpopulations
sel = c(4, 6, 7, 10, 11, 12, 13, 16)
BRO = subset(EPI, ident = sel)
saveRDS(BRO, file = "BRO.rds")


### Extraction of marker genes in the subpopulations
ALV.markers <- FindAllMarkers(ALV, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
write.table(ALV.markers, file = 'ALVMarker-genes.tsv', sep='	')

BRO.markers <- FindAllMarkers(BRO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
write.table(BRO.markers, file = 'BROMarker-genes.tsv', sep='	')



### Run monocle

ALVgene = unique(ALV.markers$gene)
BROgene = unique(BRO.markers$gene)



### ALV
# Expression Matrix 
exprs　<- GetAssayData(ALV, slot="counts",assay = "RNA")

# phenodata
pheno.data = ALV@meta.data
pheno.data = new('AnnotatedDataFrame', data = ALV@meta.data)

##feature data
genes <- data.frame(gene_short_name = rownames(ALV))
rownames(genes) <- rownames(ALV)
genes <- new('AnnotatedDataFrame', data = genes)


cdsALV <- newCellDataSet(exprs, phenoData =  pheno.data, featureData = genes, expressionFamily = uninormal())


print(dim(exprs(cdsALV)))
print(head(pData(cdsALV)))


marker_genes <- ALVgene


cdsALV <- setOrderingFilter(cdsALV, marker_genes)

cdsALV <- reduceDimension(cdsALV, norm_method = "none", 
                        reduction_method = "DDRTree",
                        max_components = 3,
                        scaling = TRUE,
                        verbose = TRUE,
						pseudo_exp = 0)

cdsALV <- orderCells(cdsALV)

saveRDS(cdsALV, file = "cdsALV.rds")




colpall = primary.colors()
colpall = sample(colpall)
colclass = c("red3", "deepskyblue3", "chocolate1")


plot_cell_trajectory(cdsALV, 
                     color_by = "Class",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colclass)


plot_cell_trajectory(cdsALV, 
                     color_by = "seurat_clusters",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colpall[1:10])

plot_cell_trajectory(cdsALV, 
                     color_by = "State",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colpall[1:9])




### BRO
# Expression Matrix 
exprs　<- GetAssayData(BRO, slot="counts",assay = "RNA")

# phenodata
pheno.data = BRO@meta.data
pheno.data = new('AnnotatedDataFrame', data = BRO@meta.data)

# feature data
genes <- data.frame(gene_short_name = rownames(BRO))
rownames(genes) <- rownames(BRO)
genes <- new('AnnotatedDataFrame', data = genes)


cdsBRO <- newCellDataSet(exprs, phenoData =  pheno.data, featureData = genes, expressionFamily = uninormal())


print(dim(exprs(cdsBRO)))
print(head(pData(cdsBRO)))


marker_genes <- BROgene


cdsBRO <- setOrderingFilter(cdsBRO, marker_genes)


cdsBRO <- reduceDimension(cdsBRO, norm_method = "none", 
                        reduction_method = "DDRTree",
                        max_components = 3,
                        scaling = TRUE,
                        verbose = TRUE,
						pseudo_exp = 0)

cdsBRO <- orderCells(cdsBRO)

saveRDS(cdsBRO, file = "cdsBRO.rds")



### Plot
plot_cell_trajectory(cdsALV, 
                     color_by = "Class",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colclass)

plot_cell_trajectory(cdsALV, 
                     color_by = "seurat_clusters",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colpall[1:10])

plot_cell_trajectory(cdsALV, 
                     color_by = "State",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colpall[1:9])



plot_cell_trajectory(cdsBRO, 
                     color_by = "Class",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5) + scale_color_manual(breaks = waiver(), values = colclass)

plot_cell_trajectory(cdsBRO, 
                     color_by = "seurat_clusters",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5) + scale_color_manual(breaks = waiver(), values = colpall[1:8])

plot_cell_trajectory(cdsBRO, 
                     color_by = "State",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1.5)  + scale_color_manual(breaks = waiver(), values = colpall[1:3])



### 3D plot
traplot = plot_cell_trajectory(cdsALV, color_by = "seurat_clusters", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 
plot_ly(x = traplot$data$data_dim_1, y = traplot$data$data_dim_2, z = traplot$data$Pseudotime, type = "scatter3d", color = traplot$data$seurat_clusters, mode = "markers", colors = colpall)

traplot = plot_cell_trajectory(cdsBRO, color_by = "seurat_clusters", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 
plot_ly(x = traplot$data$data_dim_1, y = traplot$data$data_dim_2, z = traplot$data$Pseudotime, type = "scatter3d", color = traplot$data$seurat_clusters, mode = "markers", cols = colpall)



### heatmap
diff_test_res <- differentialGeneTest(cdsALV[cdsALV@featureData@data$gene_short_name %in% ALVgene,], fullModelFormulaStr = "~Pseudotime")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
write.table(diff_test_res, "ALV-monocleGENE.tsv", sep='	')
plot_pseudotime_heatmap(cdsALV[sig_gene_names,], show_rownames = FALSE)

diff_test_res <- differentialGeneTest(cdsBRO[cdsBRO@featureData@data$gene_short_name %in% BROgene,], fullModelFormulaStr = "~Pseudotime")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
write.table(diff_test_res, "BRO-monocleGENE.tsv", sep='	')
plot_pseudotime_heatmap(cdsBRO[sig_gene_names,], show_rownames = FALSE)



### Pseudotime histgram
traALV = plot_cell_trajectory(cdsALV, color_by = "seurat_clusters", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 
traBRO = plot_cell_trajectory(cdsBRO, color_by = "seurat_clusters", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 


histALV = data.frame(cell = traALV$data$sample_name, Class = traALV$data$Class, State = traALV$data$sample_state, Seurat_cluster = traALV$data$seurat_clusters, Psuedtime = traALV$data$Pseudotime)
histBRO = data.frame(cell = traBRO$data$sample_name, Class = traBRO$data$Class, State = traBRO$data$sample_state, Seurat_cluster = traBRO$data$seurat_clusters, Psuedtime = traBRO$data$Pseudotime)

write.table(histALV, "ALV-monocleHIST.tsv", sep='	')
write.table(histBRO, "BRO-monocleHIST.tsv", sep='	')




### Movie (gif)
library(rgl)
library(magick)

traplot = plot_cell_trajectory(cdsALV, color_by = "Class", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 
colclass = c("red3", "deepskyblue3", "chocolate1")
col3d <- colclass[ as.numeric( as.factor(traplot$data$Class) ) ]

plot3d(x = traplot$data$data_dim_1, y = traplot$data$data_dim_2, z = traplot$data$Pseudotime, xlab = "Component1", ylab = "Component2", zlab = "Pseudo-time",type = "p", col = col3d)
play3d( spin3d( axis = c(0, 0, 1), rpm = 3), duration = 20 )

movie3d(
  movie = "Trajectory_ALV", 
  spin3d( axis = c(0, 0, 1), rpm = 3), duration = 20, 
  dir = "~/Analysis/COPD",
  type = "gif", 
  clean = TRUE
)



traplot = plot_cell_trajectory(cdsBRO, color_by = "Class", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 
colclass = c("red3", "deepskyblue3", "chocolate1")
col3d <- colclass[ as.numeric( as.factor(traplot$data$Class) ) ]

plot3d(x = traplot$data$data_dim_1, y = traplot$data$data_dim_2, z = traplot$data$Pseudotime, xlab = "Component1", ylab = "Component2", zlab = "Pseudo-time",type = "p", col = col3d)
play3d( spin3d( axis = c(0, 0, 1), rpm = 3), duration = 20 )

movie3d(
  movie = "Trajectory_BRO", 
  spin3d( axis = c(0, 0, 1), rpm = 3), duration = 20, 
  dir = "~/Analysis/COPD",
  type = "gif", 
  clean = TRUE
)



traplot = plot_cell_trajectory(cdsBRO, color_by = "seurat_clusters", show_branch_points = FALSE, show_tree = TRUE, cell_size = 1.5, marker_linear = TRUE) 
colpall = c("#80FF80", "#FF8080", "#FFFF00", "#8080FF", "#000000", "#008000", "#FF80FF", "#FF8000", "#000080", "#FFFF80",
			 "#FF0000", "#FF0080", "#00FF00", "#8000FF", "#80FFFF", "#FF00FF", "#808080")
col3d <- colpall[ as.numeric( as.factor(traplot$data$seurat_clusters) ) ]

plot3d(x = traplot$data$data_dim_1, y = traplot$data$data_dim_2, z = traplot$data$Pseudotime, xlab = "Component1", ylab = "Component2", zlab = "Pseudo-time",type = "p", col = col3d)
play3d( spin3d( axis = c(0, 0, 1), rpm = 3), duration = 20 )

movie3d(
  movie = "Trajectory_BRO_cluster", 
  spin3d( axis = c(0, 0, 1), rpm = 3), duration = 20, 
  dir = "~/Analysis/COPD",
  type = "gif", 
  clean = TRUE
)




