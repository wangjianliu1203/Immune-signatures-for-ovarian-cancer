
load("eSet_use.rda")
##############################################################
#PCA分析
library(openxlsx)
pData <- read.xlsx("OV_immune_nmf.xlsx", sheet = 1)

library(FactoMineR)
library(factoextra)

eSet.use <- eSet.use[,pData$Samples]
all(colnames(eSet.use)==pData$Samples)

res.pca <- PCA(X = t(eSet.use), scale.unit = F, ncp = 5, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
pdf(paste("PCA_scree_plot.pdf",sep=""),width = 6,height = 5)
scree<-fviz_eig(res.pca, addlabels = TRUE)
print(scree)
dev.off()

library(ggbiplot)
p2 <- prcomp(t(eSet.use), scale. = F)
a <- ggbiplot(p2, obs.scale = 1, var.scale = 1, groups = pData$membership2, ellipse = T, var.axes = F)+theme_light()
a <- a + ggtitle("Principal Component Analysis")
a <- a + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
a <- a + scale_color_manual(values=c("#ED0000FF", "#00468BFF"))
a <- a + ylim(c(-180, 180))
#a <- a + geom_text_repel(label = colnames(eSet.use), size = 3)
a 
ggsave(filename = "PCA.pdf", plot = a, width = 8, height = 6, dpi = 300)
dev.off()
#-------------------------------------------------------------------------------
#差异分析
library(limma)

group <- as.factor(pData$membership2)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
rownames(design) <- colnames(eSet.use)
design

contrast.matrix <- makeContrasts("Immune-non_Immune", levels = design)
contrast.matrix

fit <- lmFit(eSet.use, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

allDiff <- topTable(fit2, coef = 1, n = Inf, sort.by = "logFC")
nrDEG <- na.omit(allDiff)
head(nrDEG)
write.table(nrDEG, file="allDiff.tsv", quote=F, sep = "\t")

#输出有显著差异表达的到diffSig这个文件
diffSig <- allDiff[(abs(allDiff$logFC) > 1 & allDiff$adj.P.Val < 0.01),]
nrow(diffSig[diffSig$logFC > 0,])#230
nrow(diffSig[diffSig$logFC < 0,])#510
write.table(diffSig, file="diffSig.tsv", quote=F, sep = "\t")
#######################################################################
library(EnhancedVolcano)

q <- EnhancedVolcano(allDiff,
                     lab = rownames(allDiff),
                     x = "logFC",
                     y = "adj.P.Val",
                     ## 标记目标基因
                     selectLab = c(""),
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~-Log[10]~ 'adj-p-value'),
                     pCutoff = 0.01,
                     FCcutoff = 1,
                     pointSize = 1.0,
                     labCol = 'black',
                     labFace = 'bold',
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black',
                     ylim = c(0,28),
                     xlim = c(-3.5, 2.5),
                     colAlpha = 0.8,
                     legendLabels=c('NS',expression(Log[2]~FC),'adj-p-value',expression(adj-p-value~and~log[2]~FC)),
                     title = 'Different Expression mRNAs',
                     subtitle = paste("adj-p-value < 0.01 & |log2FC| > 1, UP DEG = ", 
                                      nrow(diffSig[diffSig$logFC > 0,]), ", DOWN DEG = ", 
                                      nrow(diffSig[diffSig$logFC < 0,]), sep = ""),
                     axisLabSize = 14, titleLabSize = 16
)
q
ggsave("DEG_EnhancedVolcano.pdf", q, width = 9, height = 9, dpi = 300)
dev.off()
########################################################################
#将差异mRNA作聚类热图
heat.mRNA <- eSet.use[rownames(eSet.use) %in% rownames(diffSig),]

library(pheatmap)
annotation <- data.frame(colnames(heat.mRNA), pData$membership2, pData$group)
row.names <- annotation[,1]
annotation <- annotation[,-1]
annotation <- as.data.frame(annotation)
rownames(annotation) <- row.names
colnames(annotation) <- c("membership", "group")
annotation$membership <- as.factor(annotation$membership)
annotation$group <- as.factor(annotation$group)
annotation_col <- list(
  group = c("1"="#00468BFF","2"="#ED0000FF","3"="#42B540FF","4"="#0099B4FF"),
  membership = c("non_Immune"="grey", "Immune"="purple")
                       )

heat.mRNA <- as.matrix(heat.mRNA)
pdf("heatmap_mRNA.pdf", onefile = F, width = 7, height = 8)
pheatmap(heat.mRNA, 
         scale = "row",
         cluster_cols = F,
         cutree_rows = 2, 
         #cutree_cols = 2,
         annotation_col = annotation, annotation_colors = annotation_col,
         color = colorRampPalette(c("#00468BFF", "white", "#ED0000FF"))(100),
         border_color = NA,
         gaps_col = c(244),
         show_rownames = F, show_colnames = F, fontsize_col = 6,fontsize =6)
dev.off()
