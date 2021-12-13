
#绘制气泡图
library(ggplot2)

slimgo <- read.delim("KEGG.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)
slimgo <- slimgo[c(1:15),]

x <- -log10(as.numeric(slimgo$PValue))
y <- factor(slimgo$Term, levels = rev(slimgo$Term))

q <- ggplot(slimgo, aes(x,y))
q <- q + geom_point(aes(size=Count, color=-log10(as.numeric(PValue)), shape=Category))
q <- q + scale_color_gradient(low ="purple", high = "yellow")
q <- q + labs(color=expression(-Log.p.value), size="Gene Count", x="-Log PValue", y="", title="DAVID Enrichment results - KEGG")
q <- q + theme_bw() + theme(
  #去网格去背景色
  axis.line = element_line(colour = "black"), axis.text = element_text(color = "black", size = 10),
  #刻度字体大小
  legend.text = element_text(size = 10), legend.title=element_text(size=10),
  #自定义气泡大小，以防太小看不见
  axis.title.x = element_text(size = 10))
q <- q + scale_size_continuous(range=c(3,7))
#q <- q + scale_shape_manual(values = c(18))
q
ggsave(filename = "enrichment_KEGG_dotplot.pdf", plot = q, device = "pdf", width = 8, height = 6, dpi = 300)
