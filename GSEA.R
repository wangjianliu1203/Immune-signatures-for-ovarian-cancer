
load("eSet_use.rda")

NMF.group <- read.delim("NMFconsensus_k2.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)
eSet.use <- eSet.use[,NMF.group$Name]

all(colnames(eSet.use)==NMF.group$Name)
write.table(eSet.use, file = "eSet_use.txt", quote = F, sep = "\t")
