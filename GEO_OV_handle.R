
#E-MTAB-2562
#读取CEL文件
library(affy)

Data <- ReadAffy()
sampleNames(Data)
N <- length(Data)

#根据cdf文件创建cdf包对数据进行标准化
library(makecdfenv)
make.cdf.package("ADXOCv1a520630.cdf", "adxocv1a520630cdf", species = "Homo sapiens", compress = TRUE)

#用RMA预处理数据
eset.rma <- rma(Data)

eSet.probe <- exprs(eset.rma)
colnames(eSet.probe) <- substr(x = colnames(eSet.probe), start = 1, stop = 9)
write.table(x = eSet.probe, file = "eSet.probe.tsv", quote = F, sep = "\t", row.names = T)
#----------------------------------------------------------------------
#对芯片进行注释
anno <- read.delim("GPL20967-3976.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)

anno_use <- anno[,c(2,3)]

eSet.probe <- as.data.frame(eSet.probe)
eSet.probe$ID_REF <- rownames(eSet.probe) 
eSet_use <- merge(eSet.probe, anno_use, by.x = "ID_REF", by.y = "ProbeSetID")
eSet_use <- eSet_use[,-1]
eSet_use <- eSet_use[nchar(eSet_use$GeneSymbol) > 0,]
eSet_use <- eSet_use[!grepl(pattern = "///", x = eSet_use$GeneSymbol),]
eSet_use <- eSet_use[!grepl(pattern = "---", x = eSet_use$GeneSymbol),]
eSet_use <- aggregate(x = eSet_use[,-ncol(eSet_use)], by = list(eSet_use[,ncol(eSet_use)]), FUN = max)
rownames(eSet_use) <- eSet_use$Group.1
eSet_use[is.na(eSet_use)] <- 0
eSet_use <- eSet_use[,-1]

write.table(x = eSet_use, file = "eSet_symbol.tsv", quote = F, sep = "\t", row.names = T)

save(eSet_use, file = "eSet_use.rda")
