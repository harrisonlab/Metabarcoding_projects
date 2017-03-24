library(data.table)
load_all("../..//metabarcoding_pipeline/scripts/myfunctions")

myfiltbiom[[1]]@sam_data$gap <- 0
myfiltbiom[[2]]@sam_data$gap <- 0

#  Pearson Correlogram
cutoff <- 17
### For cider samples - c_fix corrects a slight issue
c_fix <- function(p,o,pn){
  t1 <- plotCorrelog(p,prune_samples(sample_data(o)$block!=3,o),pn,na.add=c(9),returnCD=T)
  t2 <- plotCorrelog(p,prune_samples(sample_data(o)$block==3,o),pn,returnCD=T)
  t3 <- sapply(1:nrow(t1),function(i) 
    if(i<=nrow(t2)){cbind(rbind(t1[[i,1]],t2[[i,1]]),rbind(t1[[i,2]],t2[[i,2]]))}else{cbind(t1[[i,1]],t1[[i,2]])}
  )
  d <- as.data.frame(t(sapply(1:length(t3),function(i)
    diag(cor(t3[[i]],use="pairwise.complete.obs")[c(1,3),c(2,4)])))
  )
  d$V3 <- as.numeric(t1[[3]])
  return(d)  
}
                              
g1<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC1",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("blue","orange"),legend=F)
g2<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]],"PC1"),
		 cutoff,pc="PC1",ylim=c(-1,1),cols=c("blue","orange"),legend=T,lpos=c(0.3,0.2))
g3<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC2",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("blue","orange"),legend=F)
g4<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]],"PC2"),cutoff,pc="PC2",ylim=c(-1,1),cols=c("blue","orange"),legend=F)

pdf("correlog.pdf",width=6,height=6)
grid.arrange(
	g1+geom_text(aes(label = "A", x = 15, y = 1), color="black",size=3),
	g3+geom_text(aes(label = "C",  x = 15, y = 1),color="black",size=3),
	g2+geom_text(aes(label = "B",  x = 15, y = 1),color="black",size=3),
	g4+geom_text(aes(label = "D",  x = 15, y = 1),color="black",size=3)
)
dev.off()

g27<-plotCorrelog(data=c_fix(ITS_mypca[[2]],ITS_myfiltbiom[[2]],"PC1"),
		 cutoff,pc="PC1",ylim=c(-1,1),cols=c("blue","orange"),legend=T,lpos=c(0.3,0.3))			      
			      
lay=cbind(c(1,3),c(2,4),c(5,7),c(6,8))			      
pdf("all_correlog.pdf",width=8,height=5)
grid.arrange(
	gall[[1]]+geom_text(aes(label = "A", x = 15, y = 1), color="black",size=3),
	gall[[3]]+geom_text(aes(label = "C",  x = 15, y = 1),color="black",size=3),
	g27+geom_text(aes(label = "B",  x = 15, y = 1),color="black",size=3),
	gall[[4]]+geom_text(aes(label = "D",  x = 15, y = 1),color="black",size=3),
	g1+geom_text(aes(label = "E", x = 15, y = 1), color="black",size=3),
	g3+geom_text(aes(label = "G",  x = 15, y = 1),color="black",size=3),
	g2+geom_text(aes(label = "F",  x = 15, y = 1),color="black",size=3),
	g4+geom_text(aes(label = "H",  x = 15, y = 1),color="black",size=3),
	layout_matrix=lay
)
dev.off()			      
			      
#gt <- gtable(widths = unit(c(1, 2,3,4), "null"), heights = unit(c(1, 2, 3,4), "null"))
#gt <- gtable_add_grob(gt, ga, t = 1, b = 4, l = 4, r = 1)
#gt <- gtable_add_grob(gt, mylegend, t = 4, r = 4,l=1)
#grid.draw(gt)
