library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

##gene expression list and control
glist <- c()
universe <- c()

go <- enrichGO(glist,  OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe=universe)

reactome <- enrichPathway(bitr(glist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID)


##compare up and down regulated sets						   
glist <- list(UP=upregulated,
			DOWN = downregulated)	

go <- compareCluster(glist, fun="enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL",  universe=universe)



