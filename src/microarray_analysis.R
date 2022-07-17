setwd("W:/code/bio_project/")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)

###Fetching Data
series = "GSE48558"
platform = "GPL6244"
gset = getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
groups = c(rep("test", 13), rep("X", 27), "normal", rep("X", 3), "normal", rep("X", 23), "normal", "X", "normal", rep("X", 3), "normal", "X", rep("normal", 4), "X", "normal", rep("X", 2), rep("normal", 2), rep("X", 2), rep("normal", 2), "X", "normal", "X", "normal", "X", "normal", "X", "normal", "X", "normal", rep("X", 3), "normal",rep("X", 3), "normal", rep("X", 29), rep("normal", 7), rep("test", 2), "normal", rep("test", 3), rep("normal", 20))
sel <- which(groups != "X")
groups = groups[sel]
gset_groups = gset[, sel]
ex = exprs(gset_groups)

### Quality Control
###Box plot
pdf("results/boxplot.pdf", width = 67)
boxplot(ex)
dev.off()

###Correlation Heat map
pdf("results/cor_heatmap.pdf", width = 20, height = 20)
pheatmap(cor(ex), labels_row = groups, labels_col = groups)
dev.off()

###Principal Components Analysis
###PC on genes
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

###PC on samples
pcr = data.frame(pc$rotation[, 1:3], Group = groups)
pdf("results/pca_samples.pdf")
ggplot(pcr, aes(x=PC1, y=PC2, color = Group)) + geom_point(size= 3) + theme_bw()
dev.off()


###Differential Expression Analysis
gs <- factor(groups)
gset_groups$group <- gs
design <- model.matrix(~group + 0, gset_groups)
colnames(design) <- levels(gs)

fit <- lmFit(gset_groups, design)  # fit linear model
cts <- "test-normal"
cont.matrix <- makeContrasts(contrasts= cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="bonferroni", sort.by="logFC", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID", "adj.P.Val","logFC"))
write.table(tT, "results/AML_Normal.csv", row.names=F, sep=",", quote=FALSE)
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file="results/aml_nomral_up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file="results/aml_normal_down.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


###Correlation Heat Map for Normal Samples 
normal_groups <- c(rep("test", 13), rep("X", 27), "normal", rep("X", 3), "normal", rep("X", 23), "normal", "X", "normal", rep("X", 3), "normal", "X", rep("normal", 4), "X", "normal", rep("X", 2), rep("normal", 2), rep("X", 2), rep("normal", 2), "X", "normal", "X", "normal", "X", "normal", "X", "normal", "X", "normal", rep("X", 3), "normal",rep("X", 3), "normal", rep("X", 29), rep("normal", 7), rep("test", 2), "normal", rep("test", 3), rep("normal", 20))
normal_sel <- which(normal_groups == "normal")
normal_groups <- normal_groups[normal_sel]
normal_gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
normal_gset <- normal_gset[[idx]]
normal_gset <- normal_gset[, normal_sel]
normal_ex <- exprs(normal_gset)
normal_sourcenames <- normal_gset$source_name_ch1
pdf("results/normal_cor_heatmap.pdf", width = 15, height = 15)
pheatmap(cor(normal_ex), labels_row = normal_sourcenames, labels_col = normal_sourcenames)
dev.off()