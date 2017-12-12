X <- sumTaxa(list(as.data.frame(counts(dds,normalized=F)),taxData,colData(dds)),taxon="rank",proportional=T)
Y <- sumTaxa(list(as.data.frame(counts(dds,normalized=T)),taxData,colData(dds)),taxon="rank",proportional=T)
X$type=1
Y$type=2
X <- rbind(X,Y)

g <- ggplot(X,aes(x=factor(1),y=all,fill=rank))
g <- g + geom_bar(width = 1,stat="identity")
g <-g +facet_grid(facets=. ~ type)
g <- g + coord_polar(theta="y")
g
dev.off()
