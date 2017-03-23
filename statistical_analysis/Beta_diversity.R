```R
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

```
