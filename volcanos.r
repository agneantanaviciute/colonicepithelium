library(ggplot2)
library(ggrepel)

pdf("volcano.pdf", width=8)
ggplot(res, aes(logFC, -log10(adj.P.Val), label=labels, fill=-log10(adj.P.Val))
       ) + geom_point(aes(size=log2(Unique.peptides)), shape=21, color="black"
                      ) + geom_text_repel(size=5) + scale_fill_viridis(
                        
                      )  + labs(y="-log10 FDR", x="log2 Fold Change", 
                                size="log2 Unique Peptides", fill="-log10 FDR") + geom_hline(yintercept = -log10(0.05), lty=2, color="red")
dev.off()



