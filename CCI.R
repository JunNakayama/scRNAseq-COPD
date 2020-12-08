library(circlize)

### CCI analysis

### ligand-receptor list
### https://baderlab.org/CellCellInteractions


CCI_list <- read_csv("CCI_list.csv")
CCIgene <- read_csv("CCIgene.csv")


Chemo = CCI_list[CCI_list$AliasA %in% unlist(CCIgene),]
uni.chemo = list(Chemo[,1], Chemo[,2])
uni.chemo = unique(unlist(uni.chemo))


ALVCCI = ALV[uni.chemo,]
BROCCI = BRO[uni.chemo,]
IMMCCI = IMM[uni.chemo,]



## Reconstruction of object
new.cluster.ids <- c("AE2-1", "AE1-1", "AE2-2", "AE2-3", "AE2-4", "AE1-2", "AE1-3", "AE(2-5(Inflammation)", "AE2-6", "AE2-7")
names(new.cluster.ids) <- levels(ALVCCI)
ALVCCI <- RenameIdents(ALVCCI, new.cluster.ids)
new.cluster.ids <- c("Ciliated-1", "Ciliated-2", "Club1(AQP5+, SCGB1A1+)", "Intermediate(Ciliated, Club)", "Club(MUC5+, MUC5AC+)", "Club(SCCB1A1 high, MUC5-, MUC5AC-)", "Ciliated-3", "Basal")
names(new.cluster.ids) <- levels(BROCCI)
BROCCI <- RenameIdents(BROCCI, new.cluster.ids)
new.cluster.ids <- c("T-Cell", "NK", "CD8T", "Neutrophil-1", "Macrophage-1", "MastCell-1", "MT-1", "Macrophage-2", "Myeloid-DC", "Neutrophil-2", "MT-2", "NaiveB", 
						"MemoryB", "Plasmacytoid DC", "MastCell-2", "Macrophage-3")
names(new.cluster.ids) <- levels(IMMCCI)
IMMCCI <- RenameIdents(IMMCCI, new.cluster.ids)


CCImerge = merge(IMMCCI, c(ALVCCI, BROCCI))


### extraction of high expreesed genes ( > 2 )
unigenes = rownames(CCImerge)
mel = list()

for(i in 1: length(unigenes)){

	expr <- FetchData(object = CCImerge, vars = unigenes[i])

	if(any(expr > 2) >= 1){

		ge = CCImerge[, which(x = expr > 2)]
		ge = table(ge@active.ident)
		ho = melt(ge)
		ho = mutate(ho, Gene = unigenes[i])
		mel = rbind(mel, ho)

	}else{
		print(unigenes[i])
	}
	
}


### Convert to expressed cell ratio in the subpopulations
### Construct matrix
tocell = table(CCImerge@active.ident)
Ratiomel = list()

for(i in 1:length(tocell)){

	cell = names(tocell[i])
	LS = mel[mel[,1] == cell,]
	LS[,2] = LS[,2]/tocell[i]
	Ratiomel = rbind(Ratiomel, LS)

}

Ratiomel = dcast(Ratiomel, Var1 ~ Gene, value.var = "value")

write.table(Ratiomel, file = "CCI-RatioMAT.txt", sep = "	")



A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))


unique(list(paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|")))


ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]

CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]

CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)

a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]



CCItotal = list()

for(i in 1:nrow(CCI_perlist)){

	geneset = CCI_perlist[i,]
	Ligand = as.character(unlist(geneset[1]))
	Receptor = as.character(unlist(geneset[2]))

	LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
	ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

	for(k in 1:nrow(LigCCI)){

		Lig = LigCCI[k,3]
		score = lapply(ResCCI$value, function(x){x * Lig})
		Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
		Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

		CCItotal = bind_rows(CCItotal, Res)
	}

}



CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-1") & CCItotal$Ligand.Cell != c("MT-1"),]
CCItotal = CCItotal[CCItotal$Recepter.Cell != c("MT-2") & CCItotal$Ligand.Cell != c("MT-2"),]
GroupA = c(unlist(as.character(unique(ALVCCI@active.ident))), unlist(as.character(unique(BROCCI@active.ident))))
GroupB = c(unlist(as.character(unique(IMMCCI@active.ident))))
AtoB = CCItotal[CCItotal$Ligand.Cell %in% GroupA & CCItotal$Recepter.Cell %in% GroupB,]
BtoA = CCItotal[CCItotal$Ligand.Cell %in% GroupB & CCItotal$Recepter.Cell %in% GroupA,]




LigA = unlist(as.character(unique(AtoB$Ligand.Cell)))
AtoB.CCI = list()

for(i in 1:length(LigA)){

	G <- LigA[i]
	set = AtoB[AtoB$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		AtoB.CCI = rbind(AtoB.CCI, DF)

	}

}


LigB = unlist(as.character(unique(BtoA$Ligand.Cell)))
BtoA.CCI = list()

for(i in 1:length(LigB)){

	G <- LigB[i]
	set = BtoA[BtoA$Ligand.Cell == G, ]

	CELLset = unique(set$Recepter.Cell)

	for(k in 1:length(CELLset)){

		R <- CELLset[k]
		REset <- set[set$Recepter.Cell == R, ]
		CCIsum <- sum(REset$CCIscore)

		DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

		BtoA.CCI = rbind(BtoA.CCI, DF)

	}

}



write.table(AtoB, file = "CCI-interaction-score-AtoB.txt", sep = "	")
write.table(BtoA, file = "CCI-interaction-score-BtoA.txt", sep = "	")
write.table(AtoB.CCI, file = "CCI-interaction-score-AtoB-CELL.txt", sep = "	")
write.table(BtoA.CCI, file = "CCI-interaction-score-BtoA-CELL.txt", sep = "	")

