
#计算肿瘤TMB
library(BiocManager)
#install("PoisonAlien/TCGAmutations")

library(TCGAmutations)
library(openxlsx)

tmp <- as.data.frame(tcga_available())

dt <- tcga_load(study = "OV")
dt <- dt@data

dt1 <- as.data.frame(table(dt$Tumor_Sample_Barcode))

names(dt1) <- c('Barcode', 'Freq')

dt1$tmb <- dt1$Freq/38

names(dt1)
dt1$Barcode <- substr(dt1$Barcode, 1, 16)

write.xlsx(dt1, file = 'OV_TMB.xlsx')

groups <- read.xlsx("OV_immune_ntp.xlsx", sheet = 1)

tmb <- merge(dt1[,c(1,3)], groups[,c(1,8)], by.x = "Barcode", by.y = "Samples")
#删除离群值
tmb <- tmb[!tmb$Barcode %in% "TCGA-09-2044-01B",]

#箱线图
library(ggpubr)

p <- ggboxplot(tmb, x = "membership", y = "tmb", color = "membership", palette = c("purple", "grey"), add = "jitter", shape = "membership")
my_comparison <- list(c("1", "2"))
p <- p+stat_compare_means(comparisons = my_comparison, label = "p.signif")+stat_compare_means(label.y = max(tmb[,2])+1, label.x = 0.65)
p <- p+labs(y = paste("TMB", sep = ""))
p
ggsave(p, filename = paste("TMB_boxplot.pdf", sep = ""), device = "pdf", width = 4, height = 6)
###################################################################
# copy number variation analysis
cnv <- read.table("OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
cnv$Sample <-substr(cnv$Sample,start = 1,stop = 16)
comsam <- intersect(cnv$Sample, groups$Samples)
cnv <- cnv[which(cnv$Sample %in% comsam),]
write.table(as.data.frame(cnv),"OV_358_segment_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)

# create marker file for GISTIC2.0
marker <- cnv[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}
marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),"OV_358_marker_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

# using GenePattern to perform GISTIC2.0
library(maftools)
library(ggplot2)
library(ggpubr)

all.lesions <- "all_lesions.conf_90.txt"
amp.genes <- "amp_genes.conf_90.txt"
del.genes <- "del_genes.conf_90.txt"
scores.gis <- "scores.gistic"
gistic <- readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = F)

cnvSummary <- as.data.frame(gistic@cnv.summary)
rownames(cnvSummary) <- cnvSummary$Tumor_Sample_Barcode
comsam <- intersect(groups$Samples, rownames(cnvSummary))
cnvSummary <- merge(cnvSummary, groups[,c(1,8)], by.x = "Tumor_Sample_Barcode", by.y = "Samples")
colnames(cnvSummary)[5] <- "metacluster"
cnvSummary$metacluster <- as.factor(cnvSummary$metacluster)

# generate boxplot
my_comparisons <- list( c("0", "1"), 
                        c("0", "2"), 
                        c("1", "2"))

ggplot(data = cnvSummary,aes(x = metacluster, #分组列名
                             y = log10(Amp + 1), #连续变量列名
                             fill = metacluster))+ #按分组填充颜色
  scale_fill_manual(values = c("lightgrey", "#ED0000FF", "#42B540FF")) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("log10(Amplification)") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  stat_compare_means(comparisons = my_comparisons,method = "t.test") + 
  stat_compare_means(method = "aov", label.y = min(log10(cnvSummary$Amp + 1))-0.3)
ggsave("Figure 5C boxplot for amplification in metacluster of tcga.pdf",width = 4,height = 5)

ggplot(data = cnvSummary,aes(x = metacluster, #分组列名
                             y = log10(Del + 1), #连续变量列名
                             fill = metacluster))+ #按分组填充颜色
  scale_fill_manual(values = c("lightgrey", "#ED0000FF", "#42B540FF")) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("log10(Amplification)") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  stat_compare_means(comparisons = my_comparisons,method = "t.test") + 
  stat_compare_means(method = "aov", label.y = min(log10(cnvSummary$Del + 1))-0.3)
ggsave("Figure 5C boxplot for deletion in metacluster of tcga.pdf",width = 4,height = 5)
#######################################################################################
# seperate perform GISTIC2.0
seg.C0.cnv <- cnv[which(cnv$Sample %in% groups[groups$predict.label %in% "0",]$Samples),]
seg.C1.cnv <- cnv[which(cnv$Sample %in% groups[groups$predict.label %in% "1",]$Samples),]
seg.C2.cnv <- cnv[which(cnv$Sample %in% groups[groups$predict.label %in% "2",]$Samples),]

write.table(as.data.frame(seg.C0.cnv),"OV_metaclusterC0_segment_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.C1.cnv),"OV_metaclusterC1_segment_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.C2.cnv),"OV_metaclusterC2_segment_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)

# create marker file for GISTIC2.0
marker <- seg.C0.cnv[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}
marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),"OV_metaclusterC0_marker_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.C1.cnv[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}
marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),"OV_metaclusterC1_marker_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.C2.cnv[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}
marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),"OV_metaclusterC2_marker_forGISTIC2.0.txt",sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

# generate cytoband plot
## metaclusterC0

C0.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
pdf("Figure 5D GISTIC score for metaclusterC0 in tcga.pdf", width = 10,height = 8)
gisticChromPlot(gistic = C0.gistic, markBands = "all")
invisible(dev.off())

## metaclusterC1

C1.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
pdf("Figure 5D GISTIC score for metaclusterC1 in tcga.pdf", width = 10,height = 8)
gisticChromPlot(gistic = C1.gistic, markBands = "all")
invisible(dev.off())

## metaclusterC2

C2.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
pdf("Figure 5D GISTIC score for metaclusterC2 in tcga.pdf", width = 10,height = 8)
gisticChromPlot(gistic = C2.gistic, markBands = "all")
invisible(dev.off())
