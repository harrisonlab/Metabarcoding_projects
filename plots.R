## Volcano (like) Plots

library(calibrate)

with(res.merge,plot(log2FoldChange,log2(baseMean),pch=20, main="Volcano like plot", xlim=c(-2.5,2)))
with(subset(res.merge, padj<=.01 ), points(log2FoldChange, log2(baseMean), pch=20, col="red"))
with(subset(res.merge, abs(log2FoldChange)>=1), points(log2FoldChange, log2(baseMean), pch=20, col="orange"))
with(subset(res.merge, padj<=.01 & abs(log2FoldChange)>=1), points(log2FoldChange, log2(baseMean), pch=20, col="green"))
with(subset(res.merge, padj<=.005 & abs(log2FoldChange)>2), textxy(log2FoldChange, log2(baseMean), labs=genus, cex=.8))


# unvolcanoed
res.merge$group <-  "Not sig"
res.merge$group[res.merge$padj<=0.01] <- "sig"
g <- ggplot()
g <- g + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
g <- g + theme_bw()
g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
g <- g + scale_color_manual(values=c("black","red"))
g <- g + ylim(0,100)
g <- g + geom_point(data=res.merge, mapping=aes(x=2^log2FoldChange, y=baseMean, colour=group),size=2)
#g <- g + scale_colour_discrete(name=group)
g <- g + xlab("Fold Change")
g <- g + ylab("Base Mean")
g


res.merge$group <-  "Not sig"
res.merge$group[res.merge$padj<=0.01] <- "sig"
res.merge$group[abs(res.merge$log2FoldChange)>=1] <- "FC > 2"
res.merge$group[abs(res.merge$log2FoldChange)>=1&(res.merge$padj<=0.01)] <- "sig & FC > 2"
g <- ggplot()
g <- g + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
g <- g + theme_bw()
g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
g <- g + scale_color_manual(values=c("orange","black", "red", "green"))
g <- g + geom_point(data=res.merge, mapping=aes(x=log2FoldChange, y=log2(baseMean), colour=group),size=2)
g <- g + xlab("Log2 Fold Change")
g <- g + ylab("Log2 Base Mean")
g
dev.off()

## Correlograms

plot.corr <- function(X,sig=0.05) {
	d <- data.frame(X)
	colnames(d)<- c("pairs","distance","correlation","p")
	d$Significant<-"no"
	d$Significant[d$p<=sig]<-"yes"
	g <- ggplot(d, aes(x = distance, y = correlation,group=1,shape=Significant,colour=pairs))
	g <- g + theme_bw()
	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
	g <- g + scale_colour_gradient(low = "green", high = "red", space = "Lab", na.value = "grey50",  guide = "colourbar")
	g <- g + geom_line(na.rm=T)
	g <- g+ scale_shape_manual(values=c(1,16))
	g <- g + geom_point(na.rm=T,size=2.5,mapping=aes(shape=factor(Significant)))
	g
}


,legend.position="none"