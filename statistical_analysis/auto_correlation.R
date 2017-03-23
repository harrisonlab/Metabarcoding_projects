library(ape)
library(vegan)
library(ncf)
library(data.table)

myfiltbiom[[1]]@sam_data$gap <- 0
myfiltbiom[[2]]@sam_data$gap <- 0

cond <- "Tree"
pc.x <- lapply(seq(1,2),function(x) scores(mypca[[x]])[sample_data(myfiltbiom[[x]])$Sample==cond,])
col.x <-lapply(myfiltbiom,function(obj) sample_data(obj)[sample_data(obj)$Sample==cond,])
pc.dt <-lapply(seq(1,2),function(x) data.table(merge(pc.x[[x]],col.x[[x]],by="row.names")))

pc.reshape <- lapply(pc.dt,function(obj) {  
  y <- dcast(obj,Distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(obj)[grep("PC",names(obj))]))
  names(y)[grep("PC",names(y))] <- sub("_.*","",names(y)[grep("PC",names(y))])
  return(y)
})    


# prepare merger for biological replicates  
col.reshape <- lapply(col.x, function(obj) {
  y<-sample_data(obj)[sample_data(obj)$replicate=="a"]
  y<-y[order(y$meters)]
  return(y)
})

# Moran I test

distmat <- as.matrix(dist(cbind(sample_data(col.reshape)$distance, rep(0,24))))
distmat.inv <- 1/distmat
distmat.inv[is.infinite(distmat.inv)] <- 0
moran <- apply(pc.reshape,2,function(x) t(Moran.I(x,distmat.inv)))
names(moran)->temp
moran <- do.call(rbind,moran)
rownames(moran) <- temp
# this is similar to PCNM method shown below
distmat.bin <- (distmat > 0 $ distmat <=7.2)
moran.bin <- apply(pc.x,2,function(x) t(Moran.I(x,distmat.bin)))
rownames(moran.bin) <- rownames(moran)

# Moran correlogram

moran.mv  <- lapply(seq(1,10),function(y) correlog(sample_data(col.reshape)$distance,sample_data(col.reshape)$gap,pc.reshape$PC1,increment=y,quiet=T))
sapply(seq(1,10),function(x) plot.correlog(moran.mv[[x]]))
# I've knocked up a ggplot2 alternative plotting function, gets rid of the box around the plots
# second argument is the (two-tail) sig figure to colour points black
lapply(seq(1,10),function(x) plot.corr(moran.mv[[x]][c(1:3,5)],0.025))
dev.off()

#  Pearson Correlogram
cutoff <- 17
pc<-"PC1"
plotCorrelog(mypca[[1]],myfiltbiom[[1]],pc,cutoff=cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17))
dev.off()
d <- calcCorrelog(mypca[[1]],myfiltbiom[[1]],"PC1",c(9,17),c("Tree","Aisle"),1,F,F)

### For H samples - due to experimental design
t1 <- plotCorrelog(mypca,prune_samples(sample_data(myfiltbiom)$block!=3,myfiltbiom),pc,na.add=c(9),returnCD=T)
t2 <- plotCorrelog(mypca,prune_samples(sample_data(myfiltbiom)$block==3,myfiltbiom),pc,returnCD=T)
t3 <- sapply(1:nrow(t1),function(i) if(i<=nrow(t2)){cbind(rbind(t1[[i,1]],t2[[i,1]]),rbind(t1[[i,2]],t2[[i,2]]))}else{cbind(t1[[i,1]],t1[[i,2]])})
d <- as.data.frame(t(sapply(1:length(t3),function(i) diag(cor(t3[[i]],use="pairwise.complete.obs")[c(1,3),c(2,4)]))))
d$V3 <- as.numeric(t1[[3]])

plotCorrelog(data=d,cutoff=cutoff,pc=pc,ylim=c(-1,1))
dev.off()
