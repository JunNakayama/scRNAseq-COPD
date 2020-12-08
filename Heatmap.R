### heatmap

GENELIST = rownames(COPD)

Org = unique(COPD@active.ident)

hhh = NULL
ave = NULL
output = list()
output2 = list()


for (i in 1: length(GENELIST)){

	a = GENELIST[i]
	hoge = VlnPlot(COPD, features = a)
	hhh = data.frame(hoge$data)
		
		for(k in 1: length(Org)){

			b = Org[k]
			c = hhh[hhh$ident == b,]
			ave = mean(c[,1])

			output = rbind(output, ave)

		}

	output2 = cbind(output2, output)
	output = NULL

}

output2 = matrix(as.numeric(output2), ncol = length(GENELIST))

colnames(output2) = GENELIST
rownames(output2) = Org

write.table(output2, "AllGeneMean-matrix.tsv", sep = "	")

heatcol = brewer.pal(length(Org), "YlOrRd")
heatcol <- adjustcolor(heatcol, alpha = 1)
pheatmap(t(output2), color = heatcol)



# Extraction of GWAS genes
gene = read.delim("COPD_GWASgene.txt", header = FALSE)
gene = unlist(gene)

GWAS = subset(EPI, subset = feature = gene)
GWAS = EPI[unique(gene),]


Org = unique(GWAS@active.ident)
Org = as.factor(sort(Org))
Cls = unique(GWAS@meta.data$Class)

matlist = list()

for(i in 1:length(Org)){
	cl = Org[i]
	clcl = subset(GWAS, ident = cl)
	Idents(object = clcl) <- clcl@meta.data$Class
	avelist = AverageExpression(clcl)
	avelist = data.frame(avelist)

	matlist[i] = list(avelist)
}
names(matlist) = Org

sink("GWASmat.txt")
matlist
sink()



## Morpheus
https://software.broadinstitute.org/morpheus/
