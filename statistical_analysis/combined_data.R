#myfiltbiom <- list( 
#  dessert=prune_samples(sample_data(ITS_biom)[[11]]=="Dessert",ITS_biom),
#  cider=prune_samples(sample_data(ITS_biom)[[11]]=="Cider",ITS_biom)
#)
 
#myfiltbiom$dessert@sam_data$location <- as.factor(myfiltbiom$dessert@sam_data$meters)
#myfiltbiom$cider@sam_data$location <- as.factor(myfiltbiom$cider@sam_data$meters)

#mypca <- list(
#  dessert=plotPCA(myfiltbiom$dessert,design="1",ntop= nrow(myfiltbiom$dessert@otu_table),returnData=T,fitType="local",blind=T)
#  cider=plotPCA(myfiltbiom$cider,design="1",ntop= nrow(myfiltbiom$cider@otu_table),returnData=T,fitType="local",blind=T)
#)

myfiltbiom <- 

myfiltbiom$dessert@sam_data$location <- as.factor(myfiltbiom$dessert@sam_data$meters)

myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>5,myfiltbiom)


# get the sum of squares for tree/aisle, location and residual
sum_squares <- t(apply(mypca$x,2,function(x) t(summary(aov(x~condition+location,sample_data(myfiltbiom)))[[1]][2])))
colnames(sum_squares) <- c("condition","location","residual")
perVar <- sum_squares * mypca$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100 
