
library(data.table)
library(openxlsx)

OS.samples <- fread(input = "sp_donor.all_projects.overallSurvival_transfer_specimen.gz", sep = "\t")
OS.samples <- as.data.frame(OS.samples)
#----------------------------------------------------------------------------
#整理表达谱
eSet.sample <- fread(input = "exp_seq.OV-AU.tsv.gz", sep = "\t", )
eSet.sample <- as.data.frame(eSet.sample)
eSet.sample <- eSet.sample[,c("icgc_specimen_id", "gene_id", "raw_read_count")]

eSet <- matrix(nrow = 56206, ncol = 111)
colnames(eSet) <- levels(as.factor(eSet.sample$icgc_specimen_id))
rownames(eSet) <- levels(as.factor(eSet.sample$gene_id))
eSet <- as.data.frame(eSet)

for (i in 1:nrow(eSet.sample)) {
    eSet[eSet.sample[i,2], eSet.sample[i,1]] <- as.numeric(eSet.sample[i,3])
}
write.table(x = eSet, file = "eSet_Ensembl.txt", sep = "\t",  quote = F)
#对表达谱进行重注释
library(idmap1)

eSet$Ensembl_ID <- rownames(eSet)
anno <- annoGene(eSet$Ensembl_ID, "ENSEMBL")
write.table(anno, file = "eSet_annotation.txt", quote = F, sep = "\t", row.names = F)

eSet <- merge(eSet, anno[,c(1,3)], by.x = "Ensembl_ID", by.y = "ENSEMBL")
eSet <- eSet[,-1]
eSet <- aggregate(eSet[,-ncol(eSet)], by = list(eSet[,ncol(eSet)]), FUN = max)
rownames(eSet) <- eSet$Group.1
eSet <- eSet[,-1]

eSet.mRNA <- eSet[anno[anno$biotypes %in% "protein_coding",]$SYMBOL,]
#----------------------------------------------------------------------
#样本选择
OV.samples <- fread(input = "specimen.OV-AU.tsv.gz", sep = "\t")
OV.samples <- as.data.frame(OV.samples)

sample.sel <- OV.samples[OV.samples$specimen_type %in% "Primary tumour - solid tissue",]
write.xlsx(x = sample.sel, file = "sample_sel.xlsx")

sample.use <- intersect(x = sample.sel$icgc_specimen_id, y = OS.samples$xena_sample)
#----------------------------------------------------------------------
eSet.use <- eSet.mRNA[,sample.use]
write.table(x = eSet.use, file = "eSet_use.tsv", sep = "\t", quote = F, row.names = T)

eSet.log <- log2(eSet.use+1)
save(eSet.log, file = "eSet_use.rda")

OS.samples.use <- OS.samples[OS.samples$xena_sample %in% sample.use,]
write.table(x = OS.samples.use, file = "OS_use.tsv", sep = "\t", quote = F, row.names = F)
