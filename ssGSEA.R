
#使用estimate进行免疫纯度计算
library(estimate)
filterCommonGenes(input.f="eSet_use.txt", output.f="OV_immune_genes.gct", id="GeneSymbol")
estimateScore("OV_immune_genes.gct", "OV_estimate_score.gct", platform = "affymetrix")

#做各样本的肿瘤纯度图
plotPurity(scores="OV_estimate_score.gct", samples="all_samples")

scores <- read.delim("OV_estimate_score.gct", skip = 2, header = T)
tmp <- scores[,1]
scores <- t(scores[,3:ncol(scores)])
colnames(scores) <- tmp
scores <- as.data.frame(scores)

library(openxlsx)
write.xlsx(scores, file = "OV_estimate_score.xlsx", rowNames = T)
#------------------------------------------------------------------------
nmf.group <- read.xlsx("group.xlsx", sheet = 1) 

scores <- read.xlsx("OV_estimate_score.xlsx", sheet = 1)

immune.sel <- merge(nmf.group, scores, by.x = "Samples", by.y = "Samples")
write.xlsx(immune.sel, file = "OV_immune_nmf.xlsx", rowNames = T)
#------------------------------------------------------------------------
immune.sel <- read.xlsx("OV_immune_nmf.xlsx", sheet = 1)

library(pheatmap)

annotation <-immune.sel[,c(1,2)]
annotation <- as.data.frame(annotation)

row.names <- annotation[,1]
annotation <- annotation[,-1]
annotation <- as.data.frame(annotation)
rownames(annotation) <- row.names
annotation$annotation <- as.factor(annotation$annotation)
annotation_col <- list(annotation = c("1"="#00468BFF","2"="#ED0000FF","3"="#42B540FF","4"="#0099B4FF"))

plot <- as.data.frame(immune.sel[,c(1,3,4,5)])
rownames(plot) <- plot[,1]
plot <- plot[,-1]

table(as.factor(immune.sel[,"group"]))

pdf(file = "OV_immune_nmf.pdf", width = 12, height = 4)
pheatmap((t(plot)), 
         scale = "row",
         cluster_rows = F, cluster_cols = F, 
         show_colnames = F,
         #cutree_rows = 2,
         annotation_col = annotation, 
         annotation_colors = annotation_col,
         cellheight = 12, annotation_legend = F,
         gaps_col = c(74, 154, 223),
         color = colorRampPalette(c("#6AACCCFF", "white", "red"))(60))
dev.off()
