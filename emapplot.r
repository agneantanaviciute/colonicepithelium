library(clusterProfiler)
library(org.Hs.eg.db)

prot_res <- read.delim("proteomics.limma.results.txt")

prot_list <- list(UP = dplyr::filter(prot_res, logFC >0 & adj.P.Val < 0.05)$rownames.res.,
                  DOWN = dplyr::filter(prot_res, logFC <0 & adj.P.Val < 0.05)$rownames.res.)



prot_go2 <- enrichGO(prot_list$UP, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                    ont = c("BP"), qvalueCutoff = 0.05 )

prot_list_entrez <- bitr(prot_list$UP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

prot_go_reactome <- enrichPathway(prot_list_entrez,  qvalueCutoff = 0.5 )


geneList <- prot_res$logFC
names(geneList) <- rownames(prot_res)
geneList <- sort(geneList, decreasing = T)

map_g <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


map_g <- map_g[ !duplicated(map_g$SYMBOL),]
rownames(map_g) <- map_g$SYMBOL
map_g <- map_g[names(geneList), ]
map_g$FC <- geneList

geneList <- map_g$FC
names(geneList) <- map_g$ENTREZID

prot_go_gsea_reactome <- gsePathway(geneList, pvalueCutoff = 1)




png("proteomics_go_plot.png", res=300, height=2500, width=2500)
emapplot2(prot_go2, showCategory=30)  
dev.off()



emapplot2 <- function (x, showCategory = 30, color = "p.adjust", layout = "kk", 
          ...) 
{
  n <- enrichplot:::update_n(x, showCategory)
  geneSets <- DOSE:::geneInCategory(x)
  y <- as.data.frame(x)
  y <- y[1:n, ]
  if (n == 0) {
    stop("no enriched term found...")
  }
  else if (n == 1) {
    g <- graph.empty(0, directed = FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- y$Description
    V(g)$color <- "red"
    return(ggraph(g) + geom_node_point(color = "red", size = 5) + 
             geom_node_text(aes_(label = ~name)))
  }
  else {
    id <- y[, 1]
    geneSets <- geneSets[id]
    n <- nrow(y)
    w <- matrix(NA, nrow = n, ncol = n)
    colnames(w) <- rownames(w) <- y$Description
    for (i in 1:n) {
      for (j in i:n) {
        w[i, j] = enrichplot:::overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    wd <- melt(w)
    wd <- wd[wd[, 1] != wd[, 2], ]
    wd <- wd[!is.na(wd[, 3]), ]
    g <- graph.data.frame(wd[, -3], directed = FALSE)
    E(g)$width = sqrt(wd[, 3] * 5)
    g <- delete.edges(g, E(g)[wd[, 3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == 
                                                        y$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- y[idx, color]
    V(g)$color <- colVar
  }
  ggraph(g, layout = layout) + geom_edge_link(alpha = 0.8, 
  aes_(width = ~I(width)), colour = "darkgrey") + geom_node_point(aes_(fill = ~color, 
   size = ~size), color="black", shape=21) + geom_node_text(aes_(label = ~name), 
     repel = TRUE) + theme_void() + scale_fill_viridis() + 
    scale_size(range = c(3, 8)) + labs(fill="FDR", size="Protein Count")
}