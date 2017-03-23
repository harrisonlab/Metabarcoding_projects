library(ape)
library(vegan)
library(ncf)
library(data.table)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")

myfiltbiom[[1]]@sam_data$gap <- 0
myfiltbiom[[2]]@sam_data$gap <- 0

#  Pearson Correlogram
cutoff <- 17
pc<-"PC1"
plotCorrelog(mypca[[1]],myfiltbiom[[1]],pc,cutoff=cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17))
dev.off()

### For H samples - due to experimental design
h_fix <- function(p,o){
  t1 <- plotCorrelog(p,prune_samples(sample_data(o$block!=3,o)),pc,na.add=c(9),returnCD=T)
  t2 <- plotCorrelog(p,prune_samples(sample_data(o$block==3,o)),pc,returnCD=T)
  t3 <- sapply(1:nrow(t1),function(i) 
    if(i<=nrow(t2)){cbind(rbind(t1[[i,1]],t2[[i,1]]),rbind(t1[[i,2]],t2[[i,2]]))}else{cbind(t1[[i,1]],t1[[i,2]])}
  )
  d <- as.data.frame(t(sapply(1:length(t3),function(i)
    diag(cor(t3[[i]],use="pairwise.complete.obs")[c(1,3),c(2,4)])))
  )
  d$V3 <- as.numeric(t1[[3]])
  return(d)  
}  

plotCorrelog(data=h_fix(mypca[[2]],myfiltbiom[[2]]),cutoff=cutoff,pc=pc,ylim=c(-1,1))
dev.off()

