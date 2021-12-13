
#读取数据
library(data.table)

eSet <- fread("TCGA-OV.htseq_counts.tsv.gz", sep = "\t")
eSet <- as.data.frame(eSet)

#对表达谱进行重注释
library(idmap1)
library(openxlsx)

eSet$Ensembl_ID <- substr(eSet$Ensembl_ID, 1 ,15)
anno <- annoGene(eSet$Ensembl_ID, "ENSEMBL")
write.table(anno, file = "eSet_annotation.txt", quote = F, sep = "\t", row.names = F)

eSet <- merge(eSet, anno[,c(1,3)], by.x = "Ensembl_ID", by.y = "ENSEMBL")
eSet <- eSet[,-1]
eSet <- aggregate(eSet[,-ncol(eSet)], by = list(eSet[,ncol(eSet)]), FUN = max)
rownames(eSet) <- eSet$Group.1
eSet <- eSet[,-1]
write.table(eSet, file = "eSet_log2count_symbol.txt", quote = F, sep = "\t", row.names = T)

#提取mRNA的表达进行分析
eSet.mRNA <- eSet[anno[anno$biotypes %in% "protein_coding",]$SYMBOL,]
#---------------------------------------------------------
pData <- fread("TCGA-OV.GDC_phenotype.tsv.gz", sep = "\t")
surData <- fread("TCGA-OV.survival.tsv.gz", sep = "\t")

#将生存状态为存活0，且生存天数小于30天的样本看作随访失败
surData <- surData[!(surData$OS.time < 30 & surData$OS == 0),]

sample.use <- intersect(colnames(eSet.mRNA), pData$submitter_id.samples)
sample.use <- intersect(sample.use, surData$sample)
sample.use <- sample.use[!grepl(pattern = "-11", x = sample.use)]

eSet.use <- eSet.mRNA[,sample.use]
write.table(eSet.use, file = "eSet_use.txt", quote = F, sep = "\t", row.names = T)

pData.use <- pData[pData$submitter_id.samples %in% sample.use,]
write.xlsx(pData.use, file = "pData_use.xlsx", rowNames = F)

surData.use <- surData[surData$sample %in% sample.use,]
write.xlsx(surData.use, file = "surData_use.xlsx", rowNames = F)

