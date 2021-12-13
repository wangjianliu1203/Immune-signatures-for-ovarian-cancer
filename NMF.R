
load(file = "eSet_use.rda")

library(NMF)
library(openxlsx)
#--------------------------------------------------------------
#进行非负矩阵分解
#选择合适的rank
#可以用三个指标：cophenetic, dispersion 和silhouette
#判断最佳rank值的准则就是，cophenetic值随K变化的最大变动的前点

estim.r <- nmf(x = as.matrix(eSet.use[rowSums(eSet.use) > 0,]), rank = 2:15, nrun=10, seed=123, method = "brunet")

pdf(file = "estim.r.pdf", width = 10, height = 8)
plot(estim.r)
dev.off()

pdf(file = "cophenetic.pdf", width = 6, height = 6)
plot(2:10, estim.r$measures$cophenetic, type="b", col="purple")
dev.off()

coph <- estim.r$measures$cophenetic

coph_diff <- NULL
for (i in 2:length(coph)){
  coph_diff <- c(coph_diff, coph[i-1]-coph[i])
}
k.best <- which.max(coph_diff)+1
#-------------------------------------------------------------
#再次NMF
res <- nmf(x = as.matrix(eSet.use[rowSums(eSet.use) > 0,]), rank = k.best, nrun=50, seed = 1234, method = "brunet")

pdf(file = "basismap.pdf", width = 10, height = 10)
basismap(res)
dev.off()

group <- predict(res)
write.xlsx(x = as.data.frame(group), file = "group.xlsx", rowNames = T)

save.image("NMF.RData")

pdf(file = "consensusmap.pdf", width = 16, height = 15)
consensusmap(res, labRow = NA, labCol = NA)
dev.off()


coef(res)
basis(res)

index <- extractFeatures(res, 150L)
index
as.matrix(rownames(eSet.use)[index[[2]]])
write.xlsx(as.matrix(rownames(eSet.use)[index[[2]]]), file = "top150.xlsx")
