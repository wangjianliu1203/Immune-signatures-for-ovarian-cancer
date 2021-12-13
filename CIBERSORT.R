
library(openxlsx)

cibersort <- read.delim("CIBERSORT.Output_Job18.txt", sep = "\t", header = T, check.names = F, stringsAsFactors = F)

pData <- read.xlsx("OV_immune_ntp.xlsx", sheet = 1)

immune <- merge(cibersort, pData[,c("Samples", "predict.label")], by.x = "Input Sample", by.y = "Samples")

plot.info <- immune[,c(1,2,ncol(immune))]
plot.info$CellType <- colnames(immune)[2]
colnames(plot.info)[2] <- "Composition"

for (i in 3:(ncol(immune)-1)) {
  tmp <- immune[,c(1,i,ncol(immune))]
  tmp$CellType <- colnames(immune)[i]
  colnames(tmp)[2] <- "Composition"
  plot.info <- rbind(plot.info, tmp)
}

colnames(plot.info)[3] <- "Group"
plot.info[,3] <- as.factor(plot.info[,3])

library(ggpubr)
p <- ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition of TCGA"
) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
p
ggsave(filename = "TME.pdf", plot = p, device = "pdf", width = 10, height = 8)
dev.off()
#-------------------------------------------------------------------
q <- ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Group",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition group by Immune Group of TCGA OV",
  palette = c("lightgrey", "#ED0000FF", "#42B540FF"),
  outlier.shape = 20
) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
q
ggsave(filename = "TME_group.pdf", plot = q, device = "pdf", width = 15, height = 8)
dev.off()
#---------------------------------------------------------------------
plot.info.sel <- plot.info[plot.info$Group %in% "2",]

p2 <- ggboxplot(
  plot.info.sel,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition of TCGA"
) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
p2
ggsave(filename = "TME-group2.pdf", plot = p2, device = "pdf", width = 10, height = 8)
dev.off()
#---------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

genesel <- "Macrophages M0"

my_comparisons <- list(c("0", "1"),
                       c("0", "2"),
                       c("1", "2"))

gene.sel <- plot.info[plot.info$CellType %in% genesel,]

p <- ggplot(data = gene.sel, aes(x = Group,
                                 y = Composition,
                                 fill = Group))
p <- p + geom_boxplot(notch = T, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)

p <- p + scale_fill_manual(values = c("lightgrey", "#ED0000FF", "#42B540FF")) 
p <- p + geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                     size = 0.8, color="black") 
p <- p + geom_point(shape = 20, size=2, # 点的形状和大小
                    position = position_jitterdodge(), # 让点散开
                    color="black", alpha = 1) 
p <- p + theme_classic() + 
  ylab(genesel) +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
p <- p + stat_compare_means(method = "aov", label.y = min(plot.info$Composition)-0.2)
p
ggsave(plot = p, filename = paste0(genesel, ".pdf"), width = 6, height = 6)
#------------------------------------------------------------------------------
#使用热图进行绘制
immune.heatmap <- merge(cibersort, pData[,c("Samples", "membership","predict.label", "group")], by.x = "Input Sample", by.y = "Samples")
write.xlsx(immune.heatmap, file = "immune_heatmap.xlsx")

immune.heatmap <- read.xlsx("immune_heatmap.xlsx", sheet = 1)

library(pheatmap)

annotation <- immune.heatmap[,c(1,24:26)]
annotation <- as.data.frame(annotation)

row.names <- annotation[,1]
annotation <- annotation[,-1]
annotation <- as.data.frame(annotation)
rownames(annotation) <- row.names
annotation$membership <- as.factor(annotation$membership)
annotation$predict.label <- as.factor(annotation$predict.label)
annotation$group <- as.factor(annotation$group)

annotation_col <- list(membership = c("1"="purple", "2"="lightgrey"),
                       predict.label = c("0"="lightgrey", "1"="#ED0000FF", "2" = "#42B540FF"),
                       group = c("1"="#00468BFF","2"="#ED0000FF","3"="#42B540FF","4"="#0099B4FF"))

plot <- as.data.frame(immune.heatmap[,c(1:23)])
rownames(plot) <- plot[,1]
plot <- plot[,-1]

table(as.factor(tmp[,"predict.label"]))
table(as.factor(tmp[,"membership"]))

pdf(file = "OV_immune.pdf", width = 12, height = 12)
pheatmap((t(plot)), 
         #scale = "row",
         cluster_rows = F, cluster_cols = F, 
         show_colnames = F,
         #cutree_rows = 2,
         annotation_col = annotation, 
         annotation_colors = annotation_col,
         cellheight = 12, annotation_legend = T,
         gaps_col = c(121, 244),
         color = colorRampPalette(c("white", "pink", "red"))(60))
dev.off()

