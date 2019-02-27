
load("ABC.Rdata")
ABC@meta.data$Type[ABC@meta.data$Type == "Normal"] <- " Normal"

for ( cl in unique(ABC@ident)){
  
  if ( cl != "SKIP"){
    
    cl <-  unique(ABC@ident)[4]
    ind <-which(ABC@ident == cl)
    mast_dat <- MAST::FromMatrix(as.matrix(ABC@data[rowSums(ABC@data[, ind]) > 0, ind]), cData = ABC@meta.data[ind, ])
    zlm.output_no_pat_ct <- zlm(~ nGene  + Type, mast_dat)
    #show(zlm.output)
    
    zlm.output_pat_ct <- zlm(~ nGene  + Patient + Type, mast_dat)
    
    zlm.output_pat_ct_int2 <- zlm(~ nGene  + Type + Patient*Type + Patient, mast_dat)
    
      ``
   # show(zlm.output)
    
  }
  
  
}



ABC@meta.data$Type[ABC@meta.data$Type == "Normal"] <- " Normal"
cl <-  unique(ABC@ident)[4]
ind <-which(ABC@ident == cl)
mast_dat <- MAST::FromMatrix(as.matrix(ABC@data[rowSums(ABC@data[, ind]) > 0, ind]), cData = ABC@meta.data[ind, ])

zlm.output_no_pat_ct <- zlm(~ nGene  + Type, mast_dat)
zlm.output_pat_ct <- zlm(~ nGene  + Type + Patient, mast_dat)


summaryCond4 <- summary(zlm.output_pat_ct, doLRT='TypeInflamed') 

df <- cbind( as.data.frame(zlm.output_pat_ct_int@coefC), as.data.frame(zlm.output_no_pat_ct@coefC))
df <- df[apply(df, MARGIN=1, function(x){sum(is.na(x)) == 0}), ]

labs <- rownames(df)[which(abs(df[, 4]- df[, 12]) > .8)]

labs <- c("BLOC1S1.1", "ATP10B", "ERO1L", "PRSS3", "CRIP1", "PPP1R14A")
labs <- "MT-RNR2"
df$label <- rownames(df)
df$label[!(df$label %in% labs)] <- ""

p1 <- ggplot(df[, c(3, 4, 15)], aes(`TypeNon-inflamed`, 
                                `TypeInflamed`, fill=label, size=label)) + geom_point( shape=21, color="black"
                                  
                                ) + xlim(-3, 3) + ylim(-3, 3
                  ) + labs(title="With Patient Covariate", x="Non-Inflamed vs Healthy Coefficient", y="Inflamed vs Healthy Coefficient"
      )  + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + guides(fill=FALSE, size=FALSE)+ geom_text(aes(label=label), nudge_y = 0.5) + scale_size_discrete(range = c(1, 6)) + theme_light(base_size = 20)


p2 <- ggplot(df[, c(11, 12, 15)], aes(`TypeNon-inflamed`, 
                                  `TypeInflamed`,  fill=label, size=label)) + geom_point(shape=21, color="black"
                                    
                                  )+ xlim(-3, 3) + ylim(-3, 3)+ labs(title="Without Patient Covariate", x="Non-Inflamed vs Healthy Coefficient", y="Inflamed vs Healthy Coefficient"
                  )+ geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + guides(fill=FALSE, size=FALSE) + geom_text(aes(label=label), nudge_y = 0.5)+ scale_size_discrete(range = c(1, 6))+ theme_light(base_size = 20)



df2 <- as.data.frame(zlm.output_pat_ct_int2@coefC)
df2 <- df2[apply(df2, MARGIN=1, function(x){sum(is.na(x)) == 0}), ]

labs <- "MT-RNR2"
df2$label <- rownames(df2)
df2$label[!(df2$label %in% labs)] <- ""

p3 <- ggplot(df2[, ], aes(`TypeNon-inflamed`, 
                                      `TypeInflamed`,  fill=label, size=label)) + geom_point(shape=21, color="black"
                                                                                             
                                      )+ xlim(-3, 3) + ylim(-3, 3)+ labs(title="Without Patient Covariate", x="Non-Inflamed vs Healthy Coefficient", y="Inflamed vs Healthy Coefficient"
                                      )+ geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + guides(fill=FALSE, size=FALSE) + geom_text(aes(label=label), nudge_y = 0.5)+ scale_size_discrete(range = c(1, 6))+ theme_light(base_size = 20)





png("MT-RNR2.png", width=1300, height=1000)
cowplot::plot_grid(p1, p2, p3)
dev.off()


ggplot(data.frame(gene=ABC@data["MT-ND2", ], ABC@meta.data), aes(Patient, gene, fill=Type)) + geom_boxplot()



