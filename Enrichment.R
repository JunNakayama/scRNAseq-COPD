library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)


### Enrichment analysis of Epithelial subpopulations
h.all.genes = rownames(EPI)
h.all.genes.entrez = bitr(h.all.genes, fromType="SYMBOL", 
                         toType="ENTREZID", OrgDb="org.Hs.eg.db")
h.all.genes.entrez = h.all.genes.entrez[,2]


pmarker.entrez =  bitr(EPI.markers$gene, fromType="SYMBOL", 
                         toType="ENTREZID", OrgDb="org.Hs.eg.db")
clus = data.frame(EPI.markers$gene, EPI.markers$cluster)
pmarkers.list = merge(pmarker.entrez, clus, by = 1)
sortlist = order(pmarkers.list[,3])
pmarkers.list = pmarkers.list[sortlist,]



### GO analysis
n = unique(EPI.markers$cluster)
n = length(n)

ent.list = NULL
GORes = list()

### First, run about '0'subppulation
ent.list = pmarkers.list[pmarkers.list$EPI.markers.cluster == 0, ]
ent.list = ent.list[,2]

GOenrich <- enrichGO(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T, OrgDb = "org.Hs.eg.db")
GORes[1] = GOenrich

for(i in 1:n-1){
	ent.list = pmarkers.list[pmarkers.list$EPI.markers.cluster == i, ]
	ent.list = ent.list[,2]

	GOenrich <- enrichGO(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T, OrgDb = "org.Hs.eg.db")

	GORes[i+1] = GOenrich

}

names(GORes) = paste(unique(EPI.markers$cluster))

as.data.frame(GORes[1])


saveRDS(GORes, file = "GORes.rds")






### Pathway analysis using Reactome DB
n = unique(EPI.markers$cluster)
n = length(n)

ent.list = NULL
PathRes = list()

### First, run about '0'subppulation
ent.list = pmarkers.list[pmarkers.list$EPI.markers.cluster == 0, ]
ent.list = ent.list[,2]

GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)

PathRes[1] = GOenrich



for(i in 1:n-1){
	ent.list = pmarkers.list[pmarkers.list$EPI.markers.cluster == i, ]
	ent.list = ent.list[,2]

	GOenrich <- enrichPathway(gene = ent.list, universe = h.all.genes.entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = T)
	PathRes[i+1] = GOenrich

}

names(PathRes) = paste(unique(EPI.markers$cluster))

as.data.frame(PathRes[1])


saveRDS(PathRes, file = "PathRes.rds")


