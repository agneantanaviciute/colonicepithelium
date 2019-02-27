BatchEntropy <- function(dataset, batch0, L=100, M=100, k=500) {
  require(RANN)  
  nbatches<-length(unique(batch0))
  
  entropy<-matrix(0,L,1)
  set.seed(0) 
  for (boot in 1:L) {
    bootsamples<-sample(1:nrow(dataset),M)
    W21<-nn2(dataset,query=dataset[bootsamples,],k)
    
    for (i in 1:length(bootsamples)){
      
      for (j in 1:nbatches) {
        xi<-max(1,sum(batch0[W21$nn.idx[i,]]==j))
        entropy[boot]<-entropy[boot]+xi*log(xi)
      }
    }
  }
  
  return( (-1)*entropy/length(bootsamples) )
}


sample_size_matched <- function(seur){
  
  batches <- seur@meta.data$Sample
  bs <- c()
  for ( n in 1:length(unique(batches))){
    
    free_idents <- names(seur@ident)
    free_idents <- free_idents[!(free_idents %in% names(bs))]
    
    sample <- sample(free_idents, size = table(batches)[n], replace = F)
    
    bt <- rep(x=n, table(batches)[n])
    names(bt) <- sample
    
    bs <- c(bs, bt)
    
  }
  
  bs[names(seur@ident)]
  
}

ABC_HC@meta.data$bs1 <- sample_size_matched(ABC_HC)



batch_entropy_ABC_HC_cluster <- BatchEntropy(ABC_HC@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC_HC@ident)), M=500, k =500)
batch_entropy_ABC_HC <- BatchEntropy(ABC_HC@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC_HC@meta.data$Sample)), M=500, k =500)
batch_entropy_ABC_HC_rand <- BatchEntropy(ABC_HC@dr$tsne@cell.embeddings, batch0 = sample_size_matched(ABC_HC), M=500, k =500)


pdf("Healthy_batch_entropy.pdf")
ggplot() + geom_boxplot(aes("Positive Control\n (Clusters As Batches)", batch_entropy_ABC_HC_cluster)) + geom_boxplot(aes("Sample", batch_entropy_ABC_HC ),colour="black"
)+ geom_boxplot(aes("Negative Control\n (Random Batches)", batch_entropy_ABC_HC_rand),colour="black"
                ) + labs(x="", y="Entropy of Batch Mixing", title="HC") + theme_light(base_size = 18)
dev.off()


batch_entropy_ABC_NF_cluster <- BatchEntropy(ABC_NF@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC_NF@ident)), M=500, k =500)
batch_entropy_ABC_NF <- BatchEntropy(ABC_NF@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC_NF@meta.data$Sample)), M=500, k =500)
batch_entropy_ABC_NF_rand <- BatchEntropy(ABC_NF@dr$tsne@cell.embeddings, batch0 = sample_size_matched(ABC_NF), M=500, k =500)


pdf("Non_inflamed_batch_entropy.pdf")
ggplot() + geom_boxplot(aes("Positive Control\n (Clusters As Batches)", batch_entropy_ABC_NF_cluster)) + geom_boxplot(aes("Sample", batch_entropy_ABC_NF ),colour="black"
)+ geom_boxplot(aes("Negative Control\n (Random Batches)", batch_entropy_ABC_NF_rand),colour="black"
) + labs(x="", y="Entropy of Batch Mixing", title="Non-inflamed") + theme_light(base_size = 18)
dev.off()

batch_entropy_ABC_F_cluster <- BatchEntropy(ABC_F@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC_F@ident)), M=500, k =500)
batch_entropy_ABC_F <- BatchEntropy(ABC_F@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC_F@meta.data$Sample)), M=500, k =500)
batch_entropy_ABC_F_rand <- BatchEntropy(ABC_F@dr$tsne@cell.embeddings, batch0 = sample_size_matched(ABC_F), M=500, k =500)


pdf("Inflamed_batch_entropy.pdf")
ggplot() + geom_boxplot(aes("Positive Control\n (Clusters As Batches)", batch_entropy_ABC_F_cluster)) + geom_boxplot(aes("Sample", batch_entropy_ABC_F ),colour="black"
)+ geom_boxplot(aes("Negative Control\n (Random Batches)", batch_entropy_ABC_F_rand),colour="black"
) + labs(x="", y="Entropy of Batch Mixing", title="Inflamed") + theme_light(base_size = 18)
dev.off()


batch_entropy_ABC_cluster <- BatchEntropy(ABC@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC@ident)), M=500, k =500)
batch_entropy_ABC <- BatchEntropy(ABC@dr$tsne@cell.embeddings, batch0 = as.numeric(as.factor(ABC@meta.data$Sample)), M=500, k =500)
batch_entropy_ABC_rand <- BatchEntropy(ABC@dr$tsne@cell.embeddings, batch0 = sample_size_matched(ABC), M=500, k =500)


pdf("All_batch_entropy.pdf")
ggplot() + geom_boxplot(aes("Positive Control\n (Clusters As Batches)", batch_entropy_ABC_cluster)) + geom_boxplot(aes("Sample", batch_entropy_ABC ),colour="black"
)+ geom_boxplot(aes("Negative Control\n (Random Batches)", batch_entropy_ABC_rand),colour="black"
) + labs(x="", y="Entropy of Batch Mixing", title="All Cells") + theme_light(base_size = 18)
dev.off()

