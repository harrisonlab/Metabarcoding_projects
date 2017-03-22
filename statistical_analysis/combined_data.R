

myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom@sam_data$location <- as.factor(myfiltbiom@sam_data$meters)

tempiom <- myfiltbiom
tempiom@otu_table@.Data <-  assay(varianceStabilizingTransformation(phylo_to_des(tempiom)))

mypca <- list(
  dessert=plotPCA(tempiom,returnData=T,trans=F,filterFun=function(o,f){prune_samples(sample_data(o)[[11]]=="Dessert"&sample_data(o)[[1]]!="C",o)}),
  cider=plotPCA(tempiom,returnData=T,trans=F,filterFun=function(o,f){prune_samples(sample_data(o)[[11]]=="Cider"&sample_data(o)[[1]]!="C",o)})
)
rm(tempiom)

myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)

myfiltbiom <- list( 
  dessert=prune_samples(sample_data(myfiltbiom)[[11]]=="Dessert",myfiltbiom),
  cider=prune_samples(sample_data(myfiltbiom)[[11]]=="Cider",myfiltbiom)
)

# get the sum of squares for tree/aisle, location and residual
sum_squares <- t(apply(mypca$x,2,function(x) t(summary(aov(x~condition+location,sample_data(myfiltbiom)))[[1]][2])))
colnames(sum_squares) <- c("condition","location","residual")
perVar <- sum_squares * mypca$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100 
