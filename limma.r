library(stringr)
library(limma)

##tidy up gene names from Uniprot annotations
raw <- read.delim("raw_proteomics_data.txt")
raw$gene <- str_replace(str_split(str_match(raw[, "Description"], "GN=.+"), "\\s+", simplify = T)[, 1], "GN=", "")
rownames(raw) <- make.names(raw$gene, unique=T)

pt_dat <- raw[, c("BP1", "BP2", "BP3", "BN1", "BN2", "BN3")]
rownames(pt_dat) <- make.names(raw$gene, unique=T)

##exclude poor quality sample
pt_dat <- normalize.quantiles(as.matrix(pt_dat)[, -6])
rownames(pt_dat) <- make.names(raw$gene, unique=T)
pt_dat <- pt_dat[which(rowSums(pt_dat) > 5), ]

##make meta data table
meta_data_pt <- data.frame(BEST4=c("+", "+", "+", "-", "-"),
                           Sample=c("BP1", "BP2", "BP3", "BN1", "BN2"),
                           Donor=c("D1", "D2", "D3", "D1", "D2"))


design <- model.matrix(as.formula("~BEST4"), meta_data_pt )
fit <- limma::lmFit(log2(pt_dat + 1)[, -6], design)
ebay <- limma::eBayes(fit)
coef <- colnames(ebay$design)[2]
res <- limma::topTable(ebay, coef=coef, adjust="BH", sort.by="none", number=1e9)
res <- cbind(rownames(res), res)
res <- cbind(res, pt_dat)
table(res$adj.P.Val < 0.05)
res <- res[order(res$adj.P.Val), ]
res <- cbind(res, raw[rownames(res), c("Unique.peptides", "Anova..p.")])


res2 <- res[ res$adj.P.Val < 0.05, ]
res2 <- res2[order(res2$Unique.peptides, decreasing=T),]

