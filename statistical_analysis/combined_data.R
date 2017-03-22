library(DESeq2)
library(phyloseq)
library(ape)
library(vegan)
library(ncf)
library(data.table)
library(devtools)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")

myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom@sam_data$location <- as.factor(myfiltbiom@sam_data$meters)

tempiom <- myfiltbiom
tempiom@otu_table@.Data <-  assay(varianceStabilizingTransformation(phylo_to_des(tempiom)))

filterFun=function(o,f){
  prune_samples(sample_data(o)[[11]]==f&sample_data(o)[[1]]!="C",o)
}

mypca <- list(
  dessert=plotPCA(tempiom,returnData=T,trans=F,filterFun=filterFun,filter="Dessert"),
  cider=plotPCA(tempiom,returnData=T,trans=F,filterFun=filterFun,filter="Cider")
)
rm(tempiom)

myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)

myfiltbiom <- list( 
  dessert=filterFun(myfiltbiom,"Dessert"),
  cider=filterFun(myfiltbiom,"Cider")
)

# get the sum of squares for tree/aisle, location and residual
sum_squares <- list(
  dessert=(apply(mypca$dessert$x,2,function(x) 
    t(summary(aov(x~condition+location,data.frame(as.matrix(sample_data(myfiltbiom$dessert)))))[[1]][2]))),
  cider=(apply(mypca$cider$x,2,function(x) 
    t(summary(aov(x~condition+location,data.frame(as.matrix(sample_data(myfiltbiom$cider)))))[[1]][2])))
)
    
colnames(sum_squares) <- c("condition","location","residual")
perVar <- sum_squares * mypca$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100 
