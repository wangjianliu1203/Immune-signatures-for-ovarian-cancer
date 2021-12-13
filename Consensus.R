
library(openxlsx)

load("eSet_use.rda")

top150 <- read.xlsx("top150.xlsx", sheet = 1)

eSet.top150 <- eSet.use[top150$Top.150.genes,]
write.xlsx(eSet.top150, file = "eSet_top150.xlsx", rowNames = T)
####################################################################
immune.sel <- read.xlsx("OV_immune_nmf.xlsx", sheet = 1)
tmp <- merge(immune.sel, NMFConsensus, by.x = "Samples", by.y = "Name")
tmp <- tmp[order(tmp$membership, decreasing = F),]
write.xlsx(tmp, "OV_immune_nmf.xlsx")

tmp <- read.xlsx("OV_immune_nmf.xlsx", sheet = 1)
####################################################################
library(pheatmap)

annotation <- tmp[,c(1,2,8)]
annotation <- as.data.frame(annotation)

row.names <- annotation[,1]
annotation <- annotation[,-1]
annotation <- as.data.frame(annotation)
rownames(annotation) <- row.names
annotation$group <- as.factor(annotation$group)
annotation$membership <- as.factor(annotation$membership)
annotation_col <- list(group = c("1"="#00468BFF","2"="#ED0000FF","3"="#42B540FF","4"="#0099B4FF"),
                       membership = c("1"="purple", "2"="lightgrey"))

plot <- as.data.frame(tmp[,c(1,3,4,5)])
rownames(plot) <- plot[,1]
plot <- plot[,-1]

table(as.factor(tmp[,"membership"]))

pdf(file = "OV_immune.pdf", width = 12, height = 4)
pheatmap((t(plot)), 
         scale = "row",
         cluster_rows = F, cluster_cols = F, 
         show_colnames = F,
         #cutree_rows = 2,
         annotation_col = annotation, 
         annotation_colors = annotation_col,
         cellheight = 12, annotation_legend = F,
         gaps_col = c(244),
         color = colorRampPalette(c("#6AACCCFF", "white", "red"))(60))
dev.off()
