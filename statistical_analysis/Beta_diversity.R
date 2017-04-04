library("BiocParallel")
register(MulticoreParam(8))
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom<-prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)
uni <- UniFrac(myfiltbiom, weighted=T, normalized=T, parallel=T, fast=T)
adonis(uni~orchard*condition,as.data.frame(as.matrix(sample_data(myfiltbiom))),parallel=12,permutations=9999) 

library(GUniFrac)
library(ape)
dds <- phylo_to_des(myfiltbiom)
mytree <- phy_tree(myfiltbiom)
mytree <- root(mytree,outgroup=200,resolve.root=T)
unifracs <- GUniFrac(t(counts(dds,normalized=T)),mytree,alpha=c(0, 0.5, 1))$unifracs
dw_16 <- unifracs[, , "d_1"] # Weighted UniFrac
du_16 <- unifracs[, , "d_UW"] # Unweighted UniFrac
dv <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
d0 <- unifracs[, , "d_0"] # GUniFrac with alpha 0
d5 <- unifracs[, , "d_0.5"] # GUniFrac with alpha 0.5
adonis(as.dist(dw_16)~orchard*condition,as.data.frame(as.matrix(sample_data(myfiltbiom))),parallel=12,permutations=9999) 

dw_16 <- dw_16[row.names(sample_data(myfiltbiom)[with(sample_data(myfiltbiom),order(orchard,condition)),]),
row.names(sample_data(myfiltbiom)[with(sample_data(myfiltbiom),order(orchard,condition)),])]
plotHeatmap(dw_16)

du_16 <- du_16[row.names(sample_data(myfiltbiom)[with(sample_data(myfiltbiom),order(orchard,condition)),]),
row.names(sample_data(myfiltbiom)[with(sample_data(myfiltbiom),order(orchard,condition)),])]
plotHeatmap(du_16)

g <- plotHeatmap(dw_16)
g <-g+theme(axis.title=element_blank())
g <- g + theme(plot.margin = unit(c(1,1,2,2), "lines"))
g <- g + annotation_custom(textGrob("Grass"),xmin=35.75,xmax=35.75,ymin=-4,ymax=-4)
g <- g + annotation_custom(textGrob("Tree"),xmin=107.25,xmax=107.25,ymin=-4,ymax=-4)
g <- g + annotation_custom(textGrob("Grass"),xmin=178.75,xmax=178.75,ymin=-4,ymax=-4)
g <- g + annotation_custom(textGrob("Tree"),xmin=250.25,xmax=250.25,ymin=-4,ymax=-4)

g <- g + annotation_custom(textGrob("\nCider"),xmin=71.5,xmax=71.5,ymin=-4,ymax=-4)
g <- g + annotation_custom(textGrob("\nDessert"),xmin=214.25,xmax=214.25,ymin=-4,ymax=-4)


g <- g + annotation_custom(textGrob("\nGrass",rot=90),ymin=35.75,ymax=35.75,xmin=-12,xmax=-12)
g <- g + annotation_custom(textGrob("\nTree",rot=90),ymin=107.25,ymax=107.25,xmin=-12,xmax=-12)
g <- g + annotation_custom(textGrob("\nGrass",rot=90),ymin=178.75,ymax=178.75,xmin=-12,xmax=-12)
g <- g + annotation_custom(textGrob("\nTree",rot=90),ymin=250.25,ymax=250.25,xmin=-12,xmax=-12)

g <- g + annotation_custom(textGrob("Cider",rot=90),ymin=71.5,ymax=71.5,xmin=-12,xmax=-12)
g <- g + annotation_custom(textGrob("Dessert",rot=90),ymin=214.25,ymax=214.25,xmin=-12,xmax=-12)
gt <- ggplot_gtable(ggplot_build(g))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()

