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
colnames(sample_data(myfiltbiom))[c(1,6)] <- c("Sample","Distance")
levels(sample_data(myfiltbiom)[[1]]) <- c("C","Aisle","Grass")

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
  dessert=t(apply(mypca$dessert$x,2,function(x) 
    t(summary(aov(x~condition+location,data.frame(as.matrix(sample_data(myfiltbiom$dessert)))))[[1]][2]))),
  cider=t(apply(mypca$cider$x,2,function(x) 
    t(summary(aov(x~condition+location,data.frame(as.matrix(sample_data(myfiltbiom$cider)))))[[1]][2])))
)

sum_squares <- lapply(sum_squares,function(x) {colnames(x) <- c("condition","location","residual");x})    
    
perVar <- list(
  dessert=sum_squares$dessert * mypca$dessert$percentVar,
  cider=sum_squares$cider * mypca$cider$percentVar
)
x1 <- lapply(perVar,colSums)
x2 <- lapply(lapply(perVar,colSums),sum)
x1[[1]] <- rbind(t(data.frame(sum_sqr=x1[[1]])),"%"=x1[[1]]/x2[[1]]*100)
x1[[2]] <- rbind(t(data.frame(sum_sqr=x1[[2]])),"%"=x1[[2]]/x2[[2]]*100)

#Plot residual after removing spatial component for first couple of PCA vectors
# uncorrected
df <- lapply(mypca,function(x) {t(data.frame(t(x$x)*x$percentVar))})
# spatial removed
pc.res <- lapply(seq(1,2),function(x) resid(aov(mypca[[x]]$x~sample_data(myfiltbiom[[x]])$location)))
d <- lapply(seq(1,2),function(x) {t(data.frame(t(pc.res[[x]])*mypca[[x]]$percentVar))})
    
mygplots <- lapply(seq(1,2),function(x)    
  {g1<-plotOrd(df[[x]][,1:2],sample_data(myfiltbiom[[x]]),
          design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
          xlabel="PC1",ylabel="PC2"
  )
  g2<-plotOrd(d[[x]][,1:2],sample_data(myfiltbiom[[x]]),
          design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
          xlabel="PC1",ylabel="PC2"
  )
  g3<-plotOrd(df[[x]][,2:3],sample_data(myfiltbiom[[x]]),
          design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
          xlabel="PC2",ylabel="PC3"
  )
  g4<-plotOrd(d[[x]][,2:3],sample_data(myfiltbiom[[x]]),
          design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
          xlabel="PC2",ylabel="PC3"
  )
  return(list(g1,g2,g3,g4))
  }    
)

#Manova/Anova of first couple of pca with tree/aisle
#fit <- manova(mypca$x[,1:3]~condition,as.data.frame(as.matrix(sample_data(myfiltbiom))))
#summary(fit, test="Pillai") # could just call summary directly 
#summary(aov(mypca$x[,1]~(condition*field)+location,as.data.frame(as.matrix(sample_data(myfiltbiom)))))
#summary(aov(mypca$x[,2]~(condition*field)+location,as.data.frame(as.matrix(sample_data(myfiltbiom)))))
#summary(aov(mypca$x[,3]~(condition*field)+location,as.data.frame(as.matrix(sample_data(myfiltbiom)))))    
    
