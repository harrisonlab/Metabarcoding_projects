library(ape)
library(vegan)
library(ncf)
library(data.table)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")

myfiltbiom[[1]]@sam_data$gap <- 0
myfiltbiom[[2]]@sam_data$gap <- 0

#  Pearson Correlogram
cutoff <- 17

### For cider samples - c_fix corrects a slight issue
c_fix <- function(p,o){
  t1 <- plotCorrelog(p,prune_samples(sample_data(o)$block!=3,o),pc,na.add=c(9),returnCD=T)
  t2 <- plotCorrelog(p,prune_samples(sample_data(o)$block==3,o),pc,returnCD=T)
  t3 <- sapply(1:nrow(t1),function(i) 
    if(i<=nrow(t2)){cbind(rbind(t1[[i,1]],t2[[i,1]]),rbind(t1[[i,2]],t2[[i,2]]))}else{cbind(t1[[i,1]],t1[[i,2]])}
  )
  d <- as.data.frame(t(sapply(1:length(t3),function(i)
    diag(cor(t3[[i]],use="pairwise.complete.obs")[c(1,3),c(2,4)])))
  )
  d$V3 <- as.numeric(t1[[3]])
  return(d)  
}

pc<-"PC1"
                              
g1<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC1",cutoff=cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("blue","orange"))
g2<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]]),cutoff=cutoff,pc="PC1",ylim=c(-1,1),cols=c("blue","orange"))
g3<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC2",cutoff=cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("blue","orange"))
g4<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]]),cutoff=cutoff,pc="PC2",ylim=c(-1,1),cols=c("blue","orange"))

mylegend$layout$clip[mylegend$layout$name == "panel"] <- "off"  
  
lay=rbind(c(1,2),c(1,2),c(3,4),c(3,4),c(5,5))
 
pdf("PCA.pdf",width=6,height=6)
grid.arrange(
	g1+geom_text(aes(label = "A", x = 8, y = 6.5), color="black",size=3),
	g3+geom_text(aes(label = "C", x = 3, y = 3.5),color="black",size=3),
	g2+geom_text(aes(label = "B", x = 6.2, y = 6),color="black",size=3),
	g4+geom_text(aes(label = "D", x = 3, y = 6),color="black",size=3),
	mylegend, layout_matrix=lay
)
dev.off()
