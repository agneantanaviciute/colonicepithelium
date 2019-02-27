library(scran)
library(stringr)
library(Seurat)
##load Seurat Object
load("data.RData")

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
genes <- read.delim("genes.tsv", header=FALSE)
genes1 <- genes$V1
names(genes1) <- genes$V2
genes1 <- str_replace(genes1, "\\.\\d+", "")
cc <- cyclone( data@data, pairs=hs.pairs, verbose=T,  gene.names=genes1)
data@meta.data$phases <- cc$phases
data@meta.data$G1_score <- cc$normalized.scores$G1
data@meta.data$G2M_score <- cc$normalized.scores$G2M
data@meta.data$S_score <- cc$normalized.scores$S

##save annotated
save("data.RData")