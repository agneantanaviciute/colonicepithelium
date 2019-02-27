library("DropletUtils")

##load individual sample seurat object
load("sample_xx.RData")

sample_h5 <- list("molecule_info.h5")
s.out <- swappedDrops(sample_h5, min.frac=0.9)

iclean <- s.out$cleaned[[1]]
rownames(iclean) <- rownames(sample_XX@data)
sample_XX@raw.data <- iclean

save(sample_XX, file="sample_xx.RData")

