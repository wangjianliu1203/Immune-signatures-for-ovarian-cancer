
#瀑布图
library(maftools)

OV_maf <- read.delim("TCGA.OV.muse.somatic.maf", header = T, stringsAsFactors = F, check.names = F, sep = "\t")
OV_maf$Tumor_Sample_Barcode <- substr(OV_maf$Tumor_Sample_Barcode, 1, 16)

clinicalData <- read.delim("OV_immune_ntp.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)
OV_maf_use <- OV_maf[OV_maf$Tumor_Sample_Barcode %in% clinicalData$Tumor_Sample_Barcode,]

write.table(x = OV_maf_use, file = "TCGA.OV.somatic.maf.txt", quote = F, row.names = F, sep = "\t")
#-----------------------------------------------------------------------
OV <- read.maf(maf = "TCGA.OV.somatic.maf.txt", clinicalData = "OV_immune_ntp.txt")
getFields(OV)

#write.table(getClinicalData(OV), file = "maf_clinicalData.txt", quote = F, row.names = F, sep = "\t")
write.table(getSampleSummary(OV), file = "maf_sampleSummary.txt", quote = F, row.names = F, sep = "\t")

pdf("mafSummary.pdf", width = 10, height = 8)
plotmafSummary(maf = OV, rmOutlier = T, addStat = "median", dashboard = T, titvRaw = F)
dev.off()

#此处使用RColorBrewer的颜色，当然也可以使用任意颜色
vc_cols = RColorBrewer::brewer.pal(n = 5, name = 'Paired')
names(vc_cols) = c(
  'Nonsense_Mutation',
  'Missense_Mutation',
  'Splice_Site',
  'Translation_Start_Site',
  'Nonstop_Mutation'
)
#查看变异类型对应的颜色
print(vc_cols)

#TOP20突变的基因
fabcolors <- list(predict.label = c("0"="lightgrey", "1"="#ED0000FF", "2" = "#42B540FF"))

pdf("mafOncoplot.pdf", width = 15, height = 7)
oncoplot(maf = OV, top = 20,  colors = vc_cols, 
         clinicalFeatures = 'predict.label', sortByAnnotation = TRUE,
         annotationColor = fabcolors)
dev.off()
