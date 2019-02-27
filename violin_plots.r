library(Seurat)
library(ggplot2)
library(viridis)

##load seurat object
load("data.RData")

genes <- c("TLE4", "ELF3", "SOX4", "HES1", "MAFB")

data@ident <- factor(data@meta.data$ident)
names(data@ident) <- rownames(data@meta.data)


violin_data <- list()

for (type in unique(data@meta.data$Type)){
	df1 <- list()
	for(gene in genes){
	  
	  ind <- which(data@meta.data$Type == type)
	  
	  d <-data.frame(exp=data@data[gene, ind], cluster=data@ident[ind], gene=gene, Type=type)
	  d$cluster <- as.character(d$cluster)
	  d$cluster[which(d$cluster %in% c("Cluster 4", "Cluster 5"))] <- "Cluster 4/5" ## for goblets, merge counterpart clusters
	  
	  d$mean <- 0
	  
	  for (clust in unique(d$cluster)){
		d$mean[which(d$cluster == clust)] <- mean(d$exp[which(d$cluster == clust)])
	  }
	  
	  df1[[gene]] <- d

	}

	df <- do.call(rbind.data.frame, df1)
	violin_data[[type]] <- df
}

violin_data <- do.call(rbind, violin_data)


png("violin.png", res=300, width=4000, height=3000)
ggplot(violin_data, aes(Type, exp, fill=mean )) + geom_violin(draw_quantiles = c(0.5), scale = "width", colour="black") + labs(
  x="", y="") + theme(axis.text.y = element_blank(),
                      axis.line.y =element_blank(),
                      axis.ticks.y = element_blank(), text = element_text(size=40), axis.text.x = element_text(size=40), rect = element_blank(), 
                      strip.text.y = element_text(size=40, angle=360), strip.text.x = element_text(size=26)) + scale_fill_viridis(
                        option = "D",limits=c(0, max(violin_data$mean)),breaks=c(max(violin_data$mean)), 
                        labels=c("Avg\nExp")) + labs(fill="") + facet_grid(gene~cluster, scales="free")+ theme(axis.text.x = element_text(angle = 90, 
                                                                                                                                    hjust = 1, vjust=1)
                        ) 
dev.off()

