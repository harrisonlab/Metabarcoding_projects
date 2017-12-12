taxData<-phyloTaxaTidy(taxData,0.9,level=2)

X <- sumTaxa(list(as.data.frame(counts(dds,normalized=F)),taxData,colData(dds)),taxon="rank",proportional=T)
Y <- sumTaxa(list(as.data.frame(counts(dds,normalized=T)),taxData,colData(dds)),taxon="rank",proportional=T)
X$type=1
Y$type=2
X <- rbind(X,Y)
X$Phylum <- sub("\\(.*","",X$rank)

X <- X %>% group_by(type) %>% mutate(pos = 100-(cumsum(all)- all/2))

g <- ggplot(X,aes(x="",y=all,fill=Phylum,label = sprintf("%0.1f", round(all, digits = 1))))
g <- g + geom_bar(width = 1,stat="identity")
g <- g + scale_fill_viridis(discrete=TRUE,direction=-1)
g <- g + geom_text(data=X[X$all>4,],aes(y = pos),size=3)
#g <- g +geom_text(aes(label = sprintf("%0.2f", round(all, digits = 2))), position = position_stack(vjust = 0.5,),size=2)
g <- g + geom_label_repel(data=X[X$all<=4,],aes(y = pos), size=3, show.legend = F,colour="grey", nudge_x=0.75,nudge_y=0.75,segment.size=0.1)# + guides(fill = guide_legend(title = "rank"))
g <- g + coord_polar(theta="y")
g <- g +facet_grid(facets=. ~ type,switch="x")
g <- g + theme_blank() %+replace% theme(	
#  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  axis.text  = element_blank(),
  axis.ticks = element_blank(),
  strip.background = element_blank(),
  #strip.text.x = element_blank()
)
ggsave("FUN_pie.pdf",g+xlab(NULL)+ylab(NULL))

