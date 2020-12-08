library(network)
library(ggraph)
library(tidygraph)
library(sna)



### Normlized Closseness Centrality
### Epithelial subpopulation

cl = unique(EPI$seurat_clusters)
cl = sort(cl)
cl = as.factor(cl)
status = list()
stalist = list()
centres = list()
nomcentres = list()
sifres = list()


for(i in 1:length(cl)){

	netcell = subset(EPI, subset = seurat_clusters == cl[i])
	cNet <- cor(as.matrix(GetAssayData(netcell)))
	cNet <- 1 - cNet

	eccent = 1/apply(cNet, 2, max)
	n = nrow(cNet)
	sta = apply(cNet, 2, sum)

	status[i] = sum(sta)
	stalist[i] = list(sta)
	centres[i] = list(1/sta)
	nomcentres[i] = list((n-1)/sta)


}
names(status) = cl
names(stalist) = cl
names(centres) = cl
names(nomcentres) = cl


write.table(melt(nomcentres), file="EPIclus-normalized-closeness-centrality.tsv", sep='	')
write.table(melt(centres), file="EPIclus-closeness-centrality.tsv", sep='	')
