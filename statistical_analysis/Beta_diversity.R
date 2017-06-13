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


g <- plotHeatmap(dw_16,textSize=15)
#g_bdw, du_16,dw_16

g <- g +theme(axis.title=element_blank(),
              plot.margin = unit(c(1,1,2,2), "lines")
)

g <- g + annotation_custom(textGrob("Grass",gp = gpar(fontsize = 15)),xmin=35.75,xmax=35.75,ymin=-7,ymax=-7)
g <- g + annotation_custom(textGrob("Tree",gp = gpar(fontsize = 15)),xmin=107.25,xmax=107.25,ymin=-7,ymax=-7)
g <- g + annotation_custom(textGrob("Grass",gp = gpar(fontsize = 15)),xmin=178.75,xmax=178.75,ymin=-7,ymax=-7)
g <- g + annotation_custom(textGrob("Tree",gp = gpar(fontsize = 15)),xmin=250.25,xmax=250.25,ymin=-7,ymax=-7)
g <- g + annotation_custom(textGrob("\nCider",gp = gpar(fontsize = 15)),xmin=71.5,xmax=71.5,ymin=-9,ymax=-9)
g <- g + annotation_custom(textGrob("\nDessert",gp = gpar(fontsize = 15)),xmin=214.25,xmax=214.25,ymin=-9,ymax=-9)

g <- g + annotation_custom(textGrob("\nGrass",rot=90,gp = gpar(fontsize = 15)),ymin=35.75,ymax=35.75,xmin=-18,xmax=-18)
g <- g + annotation_custom(textGrob("\nTree",rot=90,gp = gpar(fontsize = 15)),ymin=107.25,ymax=107.25,xmin=-18,xmax=-18)
g <- g + annotation_custom(textGrob("\nGrass",rot=90,gp = gpar(fontsize = 15)),ymin=178.75,ymax=178.75,xmin=-18,xmax=-18)
g <- g + annotation_custom(textGrob("\nTree",rot=90,gp = gpar(fontsize = 15)),ymin=250.25,ymax=250.25,xmin=-18,xmax=-18)
g <- g + annotation_custom(textGrob("Cider",rot=90,gp = gpar(fontsize = 15)),ymin=71.5,ymax=71.5,xmin=-20,xmax=-20)
g <- g + annotation_custom(textGrob("Dessert",rot=90,gp = gpar(fontsize = 15)),ymin=214.25,ymax=214.25,xmin=-20,xmax=-20)

g <- g + annotation_custom(linesGrob(), xmin = 71.5, xmax = 71.5, ymin = 0, ymax = -2)
g <- g + annotation_custom(linesGrob(), xmin = 143.25, xmax = 143.25, ymin = 0, ymax = -8)
g <- g + annotation_custom(linesGrob(), xmin = 214.25, xmax = 214.25, ymin = 0, ymax = -2)

g <- g + annotation_custom(linesGrob(), ymin = 71.5, ymax = 71.5, xmin = 0, xmax = -2)
g <- g + annotation_custom(linesGrob(), ymin = 143.25, ymax = 143.25, xmin = 0, xmax = -8)
g <- g + annotation_custom(linesGrob(), ymin = 214.25, ymax = 214.25, xmin = 0, xmax = -2)

g_d <- g

g_a <- g_a + annotation_custom(textGrob("A",gp = gpar(fontsize = 19)),xmin=300,xmax=300,ymin=285,ymax=285)
g_b <- g_b + annotation_custom(textGrob("B",gp = gpar(fontsize = 19)),xmin=300,xmax=300,ymin=285,ymax=285)
g_c <- g_c + annotation_custom(textGrob("C",gp = gpar(fontsize = 19)),xmin=300,xmax=300,ymin=285,ymax=285)
g_d <- g_d + annotation_custom(textGrob("D",gp = gpar(fontsize = 19)),xmin=300,xmax=300,ymin=285,ymax=285)

gt_d <- ggplot_gtable(ggplot_build(g_d))
gt_d$layout$clip[gt_d$layout$name == "panel"] <- "off"

pdf("beta3.pdf",height=9,width=11)
grid.arrange(gt_a,gt_b,gt_c,gt_d,nrow=2,ncol=2)
#grid.draw(gt)
dev.off()

