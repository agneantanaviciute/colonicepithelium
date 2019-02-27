library(MAST)
library(Seurat)

load("data.Rdata")

##for individual clusters
cl <- "Goblets"
ind <-which(data@ident == cl)
mast_dat <- MAST::FromMatrix(as.matrix(data@data[rowSums(data@data[, ind]) > 0, ind]), cData = data@meta.data[ind, ])
    
zlm.output_no_pat_ct <- zlm(~ nGene  + Type, mast_dat)

zlm.output_pat_ct <- zlm(~ nGene  + Patient + Type, mast_dat)

summaryCond <- summary(zlm.output_pat_ct, doLRT='TypeInflamed') 

