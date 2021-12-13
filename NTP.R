
load("eSet_use.rda")

library(openxlsx)

nmf.res <- read.xlsx("OV_immune_nmf.xlsx", sheet = 1)

eSet.NTP <- eSet.use[,nmf.res[nmf.res$membership == 1,]$Samples]
write.table(eSet.NTP, file = "eSet_NTP.txt", quote = F, sep = "\t")
################################################################
ntp.res <- read.delim("NTP_sample_info.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)

tmp <- merge(nmf.res, ntp.res, by.x = "Samples", by.y = "sample.names", all.x = T)
tmp[is.na(tmp)] <- 0
write.xlsx(tmp, file = "OV_immune_ntp.xlsx")
################################################################
tmp <- read.xlsx("OV_immune_ntp.xlsx", sheet = 1)

library(pheatmap)

annotation <- tmp[,c(1,7,8)]
annotation <- as.data.frame(annotation)

row.names <- annotation[,1]
annotation <- annotation[,-1]
annotation <- as.data.frame(annotation)
rownames(annotation) <- row.names
annotation$membership <- as.factor(annotation$membership)
annotation$predict.label <- as.factor(annotation$predict.label)
annotation_col <- list(membership = c("1"="purple", "2"="lightgrey"),
                       predict.label = c("0"="lightgrey", "1"="#ED0000FF", "2" = "#42B540FF"))

plot <- as.data.frame(tmp[,c(1,2,3,4,5)])
rownames(plot) <- plot[,1]
plot <- plot[,-1]

table(as.factor(tmp[,"predict.label"]))
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
         gaps_col = c(121, 244),
         color = colorRampPalette(c("#6AACCCFF", "white", "red"))(60))
dev.off()

