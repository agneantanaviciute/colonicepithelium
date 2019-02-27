library(Seurat)

##seurat object with raw data
load("data.RData")

data <- NormalizeData(object = data,
                               normalization.method = "LogNormalize", 
                               scale.factor = 1e4)

data <- FindVariableGenes(object = data, 
                                   mean.function = ExpMean, ## defaults to mean of non-zeroes
                                   dispersion.function = LogVMR, ## default to SD
                                   x.low.cutoff = 0.0125, ## min xp
                                   x.high.cutoff = 3, ## max xp
                                   y.cutoff = 1, ## min disp
                                   do.plot=F) ## plots the dispersion/mean relationship
length(data@var.genes)



data <- ScaleData(object = data,
                           vars.to.regress = c("nUMI", "percent.mito", "Patient", "G2M_score", "G1_score"), ## regress out any variables - batches, cell cycle scores, etc.
                           do.center=T, ## 
                           do.scale=T, ##
                           model.use="linear", ## linear by default. alternatively, use glm with 'poisson' or 'negbinom' 
                           check.for.norm = T,
                           genes.use = data@var.genes) 

data <- RunPCA(object = data, 
                        pc.genes = data@var.genes, ## sustemlikeet of genes to use
                        do.print = TRUE, ## print high/low loading genes
                        pcs.print = 1:5, ## number of PCs
                        genes.print = 5,
                        pcs.compute = 50) ## number of genes for each PC


data <- ProjectPCA(object = data, 
                            do.print = T) 
PCElbowPlot(data, num.pc = 50)
data <- JackStraw(data, num.replicate = 100, do.print = T, num.pc = 50)
JackStrawPlot(object = data, PCs = 1:50)
##set PCs
nPC <- 10

data <- FindClusters(data, 
                              reduction.type = "pca", 
                              dims.use = 1:nPC, 
                              resolution = 1, 
                              print.output = 0, 
                              save.SNN = TRUE, force.recalc=TRUE)

data <- RunTSNE(object = data, 
                         dims.use = 1:nPC, 
                         do.fast = TRUE)

TSNEPlot(data,  do.return=T)

save(data, file="data.Rdata")