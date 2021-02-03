setwd("W:/code/bio_project/")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)

series = "GSE48558"
platform = "GPL6244"
gset = getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
groups = c(rep("test", 13), rep("X", 27), "normal", rep("X", 3), "normal", rep("X", 23), "normal", "X", "normal", rep("X", 3), "normal", "X", rep("normal", 4), "X", "normal", rep("X", 2), rep("normal", 2), rep("X", 2), rep("normal", 2), "X", "normal", "X", "normal", "X", "normal", "X", "normal", "X", "normal", rep("X", 3), "normal",rep("X", 3), "normal", rep("X", 29), rep("normal", 7), rep("test", 2), "normal", rep("test", 3), rep("normal", 20))
ex = exprs(gset)

pdf("results/boxplot.pdf", width = 170)
boxplot(ex)
dev.off()

###Correlation Heatmap
pdf("results/cor_heatmap.pdf", width = 25, height = 25)
pheatmap(cor(ex), labels_row = groups, labels_col = groups)
dev.off()

###Principal Components Analysis
pc = prcomp(ex)
pdf("results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

ex_scaled = t(scale(t(ex), scale = FALSE))
pc = prcomp(ex_scaled)
pdf("results/pc_scaled.pdf")
plot(pc)
plot(pc$x[, 1:2])
dev.off()

pcr = data.frame(pc$rotation[, 1:3], Group = groups)
pdf("results/pca_samples.pdf")
ggplot(pcr, aes(x=PC1, y=PC2, color = Group)) + geom_point(size= 3) + theme_bw()
dev.off()


###Differential Expression Analysis
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="bonferroni", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
