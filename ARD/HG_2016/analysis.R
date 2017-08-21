#===============================================================================
#       Load libraries
#===============================================================================

library(phyloseq)
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(dplyr)
library(plyr)
library(devtools)
load_all("../../metabarcoding_pipeline/scripts/myfunctions")
library(vegan)
library(ape)
library(ggplot2)
library(reshape2)


#===============================================================================
#       Load data 
#===============================================================================

get_biom <- function(biom,colData) {
	X<-import_biom(biom)
	sample_data(X) <- read.table(colData,header=T,sep="\t",row.names=1)
	tax_table(X) <- phyloTaxaTidy(tax_table(X),0.65)
	return(X)
}

biom16 <- get_biom("16S.taxa.biom","colData")
biomITS <- get_biom("ITS.taxa.biom","colData")

mybioms <- list(Bacteria=biom16,Fungi=biomITS)


#===============================================================================
#       Pool Data
#===============================================================================
myfiltbioms <- lapply(mybioms,function(obj) prune_samples(	
	(sample_data(obj)[[10]]!="duplicate")&
	(sample_data(obj)[[1]]!="C")&
	(colSums(otu_table(obj))>999)
,obj))


myfiltbioms <- lapply(myfiltbioms,function(obj) {
	obj@sam_data$location <- as.factor(obj@sam_data$meters)
	colnames(sample_data(obj))[c(1,6,11)] <- c("Sample","Distance","Orchard")
	levels(sample_data(obj)[[1]]) <- c("Aisle","Tree")
	return(obj)
})


Ldds <- lapply(myfiltbioms,function(obj) phylo_to_des(obj,fit=F,calcFactors=geoSet))

Ldds$Fungi$Group <- as.factor(paste(Ldds$Fungi$Orchard,Ldds$Fungi$Sample,Ldds$Fungi$location,sep="_"))
Ldds$Bacteria$Group <- as.factor(paste(Ldds$Bacteria$Orchard,Ldds$Bacteria$Sample,Ldds$Bacteria$location,sep="_"))
Ldds<-lapply(Ldds,function(o) collapseReplicates2(o,groupby=o$Group))
#sizeFactors(Ldds$Fungi)<-1
#sizeFactors(Ldds$Bacteria)<-1
Ldds$Fungi$Orchard<- droplevels(Ldds$Fungi$Orchard)
Ldds$Fungi$Sample<- droplevels(Ldds$Fungi$Sample)
Ldds$Fungi$location<- droplevels(Ldds$Fungi$location)
Ldds$Bacteria$Orchard<- droplevels(Ldds$Bacteria$Orchard)
Ldds$Bacteria$Sample<- droplevels(Ldds$Bacteria$Sample)
Ldds$Bacteria$location<- droplevels(Ldds$Bacteria$location)

levels(Ldds$Fungi$Sample)<-c("Grass","Tree")
levels(Ldds$Bacteria$Sample)<-c("Grass","Tree")

#countData <- lapply(Ldds,function(o) assay(o))
countData <- lapply(Ldds,function(o) counts(o,normalize=T))
colData <- lapply(Ldds,function(o) o@colData)

#===============================================================================
#       Base analysis
#===============================================================================

d<-rowSums(counts(Ldds[[1]][,Ldds[[1]]@colData$Orchard=="Dessert"&Ldds[[1]]@colData$Sample=="Tree"],normalized=T))
d<-cumsum(d[order(d,decreasing=T)])

dd<-rowSums(assay(Ldds[[1]][,Ldds[[1]]@colData$Orchard=="Dessert"&Ldds[[1]]@colData$Sample=="Grass"],normalized=T))
dd<-cumsum(dd[order(dd,decreasing=T)])

ddd<-rowSums(counts(Ldds[[1]][,Ldds[[1]]@colData$Orchard=="Cider"&Ldds[[1]]@colData$Sample=="Tree"],normalized=T))
ddd<-cumsum(ddd[order(ddd,decreasing=T)])

dddd<-rowSums(assay(Ldds[[1]][,Ldds[[1]]@colData$Orchard=="Cider"&Ldds[[1]]@colData$Sample=="Grass"],normalized=T))
dddd<-cumsum(dddd[order(dddd,decreasing=T)])

df <- cbind(d,dd,ddd,dddd)
colnames(df) <- c("Dessert Tree","Dessert Grass","Cider Tree", "Cider Grass")
rownames(df) <- seq(1,nrow(df))
myline<-apply(df,2,function(o) length(o)-length(o[o>=(max(o,na.rm=T)*0.8)]))
df<- melt(log10(df),id=rownames)
colnames(df)[1:2]<-c("OTUS","Location")

pdf("POOLED_OTU_counts_n2.pdf")


g <- ggplot(data=df,aes(x=OTUS,y=value,colour=Location))
g <- g + theme_bw()
g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(axis.line.x = element_line(size=0.3,colour = "black"),axis.line.y = element_line(size=0.3,colour = "black"),axis.text = element_text(colour = "black"),text=element_text(size=16))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
g+geom_line(size=1.5)+scale_colour_manual(values=cbbPalette)+ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")+geom_vline(xintercept=myline,colour=cbbPalette[1:4])

i=2
d<-rowSums(assay(Ldds[[i]][,Ldds[[i]]@colData$Orchard=="Dessert"&Ldds[[i]]@colData$Sample=="Tree"],normalized=T))
d<-cumsum(d[order(d,decreasing=T)])

dd<-rowSums(assay(Ldds[[i]][,Ldds[[i]]@colData$Orchard=="Dessert"&Ldds[[i]]@colData$Sample=="Grass"],normalized=T))
dd<-cumsum(dd[order(dd,decreasing=T)])

ddd<-rowSums(assay(Ldds[[i]][,Ldds[[i]]@colData$Orchard=="Cider"&Ldds[[i]]@colData$Sample=="Tree"],normalized=T))
ddd<-cumsum(ddd[order(ddd,decreasing=T)])

dddd<-rowSums(assay(Ldds[[i]][,Ldds[[i]]@colData$Orchard=="Cider"&Ldds[[i]]@colData$Sample=="Grass"],normalized=T))
dddd<-cumsum(dddd[order(dddd,decreasing=T)])
len<-min(sapply(list(d,dd,ddd,dddd),length))
df <- cbind(d[1:len],dd[1:len],ddd[1:len],dddd[1:len])
colnames(df) <- c("Dessert Tree","Dessert Grass","Cider Tree", "Cider Grass")
rownames(df) <- seq(1,nrow(df))
myline<-apply(df,2,function(o) length(o)-length(o[o>=(max(o,na.rm=T)*0.8)]))

df<- melt(log10(df),id=rownames)
colnames(df)[1:2]<-c("OTUS","Location")

g <- ggplot(data=df,aes(x=OTUS,y=value,colour=Location))
g <- g + theme_bw()
g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(axis.line.x = element_line(size=0.3,colour = "black"),axis.line.y = element_line(size=0.3,colour = "black"),axis.text = element_text(colour = "black"),text=element_text(size=16),legend.position="none")
g<-g+geom_line(size=1.5)+scale_colour_manual(values=cbbPalette)+ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")+geom_vline(xintercept=myline,colour=cbbPalette[1:4])+coord_cartesian(xlim = c(0, 100)) 

dev.off()


## plots


lapply(res.merge,function(obj) {
	with(obj,plot(log2FoldChange,log10(baseMean),pch=20, xlim=c(-6,6),bty="n",
	xlab=expression("Log"[2]*" Fold Change"),ylab=expression("Log"[10]*" Mean Expression")))
	with(subset(obj, padj<0.5 ), points(log2FoldChange, log10(baseMean), pch=20, col="#E69F00"))
	with(subset(obj, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="#56B4E9"))
	with(subset(obj, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="#009E73"))
})
dev.off()


## taxonomy calculations

test <- as.data.frame(as.matrix(tax_table(myfiltbioms[[1]])[,8:14]))
test<- as.matrix(as.data.frame(lapply(test,function(i) as.numeric(levels(i)[i]))))
apply(test,2,function(x) sum(x>=0.65)/length(x)*100)


#===============================================================================
#       PCA analysis
#===============================================================================

LddsFilt <- lapply(Ldds,function(o) o[rowSums(counts(o,normalize=T))>5,])

pca_des <- function(obj) {
	X<-prcomp(t(assay(varianceStabilizingTransformation(obj))))
	X$percentVar<-X$sdev^2/sum(X$sdev^2)
	X
}

mypca <- lapply(LddsFilt,pca_des)
	

#pca anova
#Bacteria
lapply(seq(1:4),function(x)
	summary(aov(mypca[[1]]$x[,x]~location+(Orchard*Sample),colData[[1]]))
)


lapply(seq(1:4),function(x)
	summary(aov(mypca[[1]]$x[,x]~location+(Orchard*Sample),colData[[1]]))[[1]][[2]]/
	sum(summary(aov(mypca[[1]]$x[,x]~location+(Orchard*Sample),colData[[1]]))[[1]][[2]])*100
)


#Fungi
lapply(seq(1:4),function(x)
	summary(aov(mypca[[2]]$x[,x]~location+(Orchard*Sample),colData[[2]]))
)

lapply(seq(1:4),function(x)
	summary(aov(mypca[[2]]$x[,x]~location+(Orchard*Sample),colData[[2]]))[[1]][[2]]/
	sum(summary(aov(mypca[[2]]$x[,x]~location+(Orchard*Sample),colData[[2]]))[[1]][[2]])*100
)

# Calculate sum of squares
#Bacteria
sum_squares <- t(apply(mypca[[1]]$x,2,function(x) 
  t(summary(aov(x~location+(Orchard*Sample),colData[[1]]))[[1]][2]))
)
colnames(sum_squares) <- c("location","orchard","condition","orchard:sample","residual")
x<-t(apply(sum_squares,1,prop.table))
perVar <- x * mypca[[1]]$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100
#Fungi
sum_squares <- t(apply(mypca[[2]]$x,2,function(x) 
  t(summary(aov(x~location+(Orchard*Sample),colData[[2]]))[[1]][2]))
)
colnames(sum_squares) <- c("location","orchard","condition","orchard:sample","residual")
x<-t(apply(sum_squares,1,prop.table))
perVar <- x * mypca[[2]]$percentVar
colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100

# PCA plots

df <- lapply(mypca,function(o) t(data.frame(t(o$x)*o$percentVar)))

g1 <- plotOrd(df[[2]],colData[[2]],shapes=c("Orchard","Sample"),design="Distance",xlabel="PC1",ylabel="PC2",continuous=T,dimx=1,dimy=2,colourScale=c("black","lightblue"),legend=F,textSize=14)
g2 <- plotOrd(df[[1]],colData[[1]],shapes=c("Orchard","Sample"),design="Distance",xlabel="PC1",ylabel="PC2",continuous=T,dimx=1,dimy=2,colourScale=c("black","lightblue"),legend=F,textSize=14)

g2 <- g2 + theme(legend.direction="horizontal", 
		 legend.position="bottom",
		 legend.justification=c(0,0),
		 legend.box="vertical",
		 legend.box.just="left",
		 legend.text=element_text(size=12),
		 legend.title=element_text(size=12),
		 axis.line.x = element_line(size=0.3,colour = "black"),
		 axis.line.y = element_line(size=0.3,colour = "black"),
		 axis.text = element_text(colour = "black"),
		 plot.margin=unit(c(2.5,1,0.5,0.5), "lines")
)	

### Main figure		  
pdf("POOLED_ITS_16S_orchards2.pdf",width=7,height=6.5)	
g3 <- ggplotGrob(g1+annotate("text",label=paste("A"), x=25, y=3,size=5))
g4<-  ggplotGrob(g2+annotate("text",label=paste("B"), x=60, y=8,size=5))
g <- rbind(g3, g4, size="first")
#grid.newpage()
grid.draw(g)
dev.off()


### PCA per orchard
#phylist <- lapply(seq(1,2),function(i) phyloseq(otu_table(counts(Ldds[[i]],normalize=F),taxa_are_rows=T),sample_data(as.data.frame(colData[[i]])),tax_table(tax_table(myfiltbiom[[i]]))))
#tempiom <- phylist[[2]]
#tempiom@otu_table@.Data <-  assay(varianceStabilizingTransformation(phylo_to_des(tempiom)))

#..then as per git (mostly)

#filterFun=function(o,f){
#  o<-prune_samples(sample_data(o)[[11]]==f,o)
#  prune_taxa(rowSums(otu_table(o))>5,o)
#}	

#phylist <- lapply(phylist,function(o) prune_taxa(rowSums(otu_table(o))>5,o))

#myfiltbiom<-phylist[[2]]

#mypca <- 	prcomp(t(otu_table(myfiltbiom)))
#mypca$percentVar<-mypca$sdev^2/sum(mypca$sdev^2)


phylist <- lapply(seq(1,2),function(i) phyloseq(
	otu_table(counts(LddsFilt[[i]],normalize=T),taxa_are_rows=T),
	sample_data(as.data.frame(colData[[i]]))
	#tax_table(myfiltbioms[[i]])
))

mypca <- list(
  Bacteria_Dessert = pca_des(LddsFilt$Bacteria[,LddsFilt$Bacteria$Orchard=="Dessert"]), 
  Bacteria_Cider =   pca_des(LddsFilt$Bacteria[,LddsFilt$Bacteria$Orchard=="Cider"]),
  Fungi_Dessert =	 pca_des(LddsFilt$Fungi[,LddsFilt$Fungi$Orchard=="Dessert"]),
  Fungi_Cider = 	 pca_des(LddsFilt$Fungi[,LddsFilt$Fungi$Orchard=="Cider"])
)


myfiltbiom <- list(
	Bacteria_Dessert=prune_samples(sample_data(phylist[[1]])$Orchard=="Dessert",phylist[[1]]),
	Bacteria_Cider=prune_samples(sample_data(phylist[[1]])$Orchard=="Cider",phylist[[1]]),
	Fungi_Dessert=prune_samples(sample_data(phylist[[2]])$Orchard=="Dessert",phylist[[2]]),
	Fungi_Cider=prune_samples(sample_data(phylist[[2]])$Orchard=="Cider",phylist[[2]])
)

df <- lapply(mypca,function(o) t(data.frame(t(o$x)*o$percentVar)))

# location effect removed...
pc.res <- lapply(seq(1,4),function(i) resid(aov(mypca[[i]]$x~myfiltbiom[[i]]@sam_data$location)))
d <- lapply(seq(1,4),function(i) t(data.frame(t(pc.res[[i]])*mypca[[i]]$percentVar)))


pdf("bd1.pdf")
g1<-plotOrd(df[[1]],sample_data(myfiltbiom[[1]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g1;dev.off()
pdf("bd2.pdf")
g2<-plotOrd(d[[1]],sample_data(myfiltbiom[[1]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g2;dev.off()
pdf("bc1.pdf")
g3<-plotOrd(df[[2]],sample_data(myfiltbiom[[2]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g3;dev.off()
pdf("bc2.pdf")
g4<-plotOrd(d[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g4;dev.off()
pdf("fd1.pdf")
g5<-plotOrd(df[[3]],sample_data(myfiltbiom[[3]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g5;dev.off()
pdf("fd2.pdf")
g6<-plotOrd(d[[3]],sample_data(myfiltbiom[[3]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g6;dev.off()
pdf("fc1.pdf")
g7<-plotOrd(df[[4]],sample_data(myfiltbiom[[4]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g7;dev.off()
pdf("fc2.pdf")
g8<-plotOrd(d[[4]],sample_data(myfiltbiom[[4]]),
	design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
	cluster=0.95,centers=1,xlabel="PC1",ylabel="PC2",legend=F,textSize=16,xlim=c(-10,10),ylim=c(-10,10)
)
g8;dev.off()



g1;g2;g3;g4;g5;g6;g7;g8
dev.off()

mylegend <- ggplot_legend( plotOrd(df[[2]][,1:2],sample_data(myfiltbiom[[2]]),
	  design="Distance",shapes="Sample",continuous=T,colourScale=c("black","lightblue"),
    )	+theme(legend.direction="horizontal",legend.box="horizontal")
)

#===============================================================================
#       autocorrelation
#===============================================================================
LddsFilt <- lapply(Ldds,function(o) o[rowSums(counts(o,normalize=T))>5,])

phylist <- lapply(seq(1,2),function(i) phyloseq(
	otu_table(counts(LddsFilt[[i]],normalize=T),taxa_are_rows=T),
	sample_data(as.data.frame(colData[[i]]))
	#tax_table(myfiltbioms[[i]])
))


mypca <- list(
  Bacteria_Dessert = pca_des(LddsFilt$Bacteria[,LddsFilt$Bacteria$Orchard=="Dessert"]), 
  Bacteria_Cider =   pca_des(LddsFilt$Bacteria[,LddsFilt$Bacteria$Orchard=="Cider"]),
  Fungi_Dessert =	 pca_des(LddsFilt$Fungi[,LddsFilt$Fungi$Orchard=="Dessert"]),
  Fungi_Cider = 	 pca_des(LddsFilt$Fungi[,LddsFilt$Fungi$Orchard=="Cider"])
)


myfiltbiom <- list(
	Bacteria_Dessert=prune_samples(sample_data(phylist[[1]])$Orchard=="Dessert",phylist[[1]]),
	Bacteria_Cider=prune_samples(sample_data(phylist[[1]])$Orchard=="Cider",phylist[[1]]),
	Fungi_Dessert=prune_samples(sample_data(phylist[[2]])$Orchard=="Dessert",phylist[[2]]),
	Fungi_Cider=prune_samples(sample_data(phylist[[2]])$Orchard=="Cider",phylist[[2]])
)


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
		       
g1<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC1",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g2<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]],"PC1"), cutoff,pc="PC1",ylim=c(-1,1),cols=c("black","lightblue"),legend=F)
g3<-plotCorrelog(mypca[[1]],myfiltbiom[[1]],"PC2",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g4<-plotCorrelog(data=c_fix(mypca[[2]],myfiltbiom[[2]],"PC2"),cutoff,pc="PC2",ylim=c(-1,1),cols=c("black","lightblue"),legend=F)

g11<-plotCorrelog(mypca[[3]],myfiltbiom[[3]],"PC1",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g12<-plotCorrelog(data=c_fix(mypca[[4]],myfiltbiom[[4]],"PC1"), cutoff,pc="PC1",ylim=c(-1,1),cols=c("black","lightblue"),legend=T,lpos=c(0.275,0.25))
g13<-plotCorrelog(mypca[[3]],myfiltbiom[[3]],"PC2",cutoff,xlim=NULL,ylim=c(-1,1),na.add=c(9,17),cols=c("black","lightblue"),legend=F)
g14<-plotCorrelog(data=c_fix(mypca[[4]],myfiltbiom[[4]],"PC2"),cutoff,pc="PC2",ylim=c(-1,1),cols=c("black","lightblue"),legend=F)			      
		      
lay=cbind(c(1,3),c(2,4),c(5,7),c(6,8))			      
pdf("all_correlog_bb2.pdf",width=8,height=5)
grid.arrange(
	g11+geom_text(aes(label = "A", x = 15, y = 1), color="black",size=3),
	g13+geom_text(aes(label = "C",  x = 15, y = 1),color="black",size=3),
	g12+geom_text(aes(label = "B",  x = 15, y = 1),color="black",size=3),
	g14+geom_text(aes(label = "D",  x = 15, y = 1),color="black",size=3),
	g1+geom_text(aes(label = "E", x = 15, y = 1), color="black",size=3),
	g3+geom_text(aes(label = "G",  x = 15, y = 1),color="black",size=3),
	g2+geom_text(aes(label = "F",  x = 15, y = 1),color="black",size=3),
	g4+geom_text(aes(label = "H",  x = 15, y = 1),color="black",size=3),
	layout_matrix=lay
)
dev.off()			      

#===============================================================================
#       diversity analysis
#===============================================================================

##### ALPHA Diversity ######

alpha_counts <- lapply(countData,function(X) {
	X[X==0]<-NA
	X[X<1] <- 1
	X[is.na(X)] <- 0
	return(round(X,0))
})


phylist <- lapply(seq(1,2),function(i) phyloseq(otu_table(alpha_counts[[i]],taxa_are_rows=T),sample_data(as.data.frame(colData[[i]])),tax_table(tax_table(myfiltbioms[[i]]))))

all_alpha <- lapply(phylist,function(o) plot_richness(o,returnData=T))


data.frame(prop.table(summary(aov(Chao1~location+(Sample*Orchard),all_alpha[[2]]))[[1]][c(2)])*100,
           summary(aov(Chao1~location+(Sample*Orchard),all_alpha[[2]]))[[1]][5])
data.frame(prop.table(summary(aov(Shannon~location+(Sample*Orchard),all_alpha[[2]]))[[1]][c(2)])*100,
           summary(aov(Shannon~location+(Sample*Orchard),all_alpha[[2]]))[[1]][5])
data.frame(prop.table(summary(aov(Simpson~location+(Sample*Orchard),all_alpha[[2]]))[[1]][c(2)])*100,
           summary(aov(Simpson~location+(Sample*Orchard),all_alpha[[2]]))[[1]][5])
           
           
data.frame(prop.table(summary(aov(Chao1~location+(Sample*Orchard),all_alpha[[1]]))[[1]][c(2)])*100,
           summary(aov(Chao1~location+(Sample*Orchard),all_alpha[[1]]))[[1]][5])
data.frame(prop.table(summary(aov(Shannon~location+(Sample*Orchard),all_alpha[[1]]))[[1]][c(2)])*100,
           summary(aov(Shannon~location+(Sample*Orchard),all_alpha[[1]]))[[1]][5])
data.frame(prop.table(summary(aov(Simpson~location+(Sample*Orchard),all_alpha[[1]]))[[1]][c(2)])*100,
           summary(aov(Simpson~location+(Sample*Orchard),all_alpha[[1]]))[[1]][5])           


summary(aov(Chao1~location+(Sample*Orchard),all_alpha[[2]]))
summary(aov(Shannon~location+(Sample*Orchard),all_alpha[[2]]))
summary(aov(Simpson~location+(Sample*Orchard),all_alpha[[2]]))
           

summary(aov(Chao1~location+(Sample*Orchard),all_alpha[[1]]))
summary(aov(Shannon~location+(Sample*Orchard),all_alpha[[1]]))
summary(aov(Simpson~location+(Sample*Orchard),all_alpha[[1]]))



#### CHANGE i ####
i=1 #or 2
sample_data(phylist[[i]])$Class <- paste(sample_data(phylist[[i]])$Orchard,sample_data(phylist[[i]])$Sample,sep=" ")
sample_data(phylist[[i]])$Class[sample_data(phylist[[i]])$Class=="Cider Grass"] <- "C-G"
sample_data(phylist[[i]])$Class[sample_data(phylist[[i]])$Class=="Cider Tree"] <- "C-T"
sample_data(phylist[[i]])$Class[sample_data(phylist[[i]])$Class=="Dessert Grass"] <- "D-G"
sample_data(phylist[[i]])$Class[sample_data(phylist[[i]])$Class=="Dessert Tree"] <- "D-T"


g1 <- plot_richness(phylist[[2]],x="Class",color="Distance",measures=c("Chao1", "Shannon", "Simpson"))
g2 <- plot_richness(phylist[[1]],x="Class",color="Distance",measures=c("Chao1", "Shannon", "Simpson"))

title.A <- textGrob(label = "A",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = 0,gp = gpar(fontsize = 16))
title.B <- textGrob(label = "B",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = -1,gp = gpar(fontsize = 16))

#dev.off()
g1 <- g1 + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 plot.margin=unit(c(-1,1,2.5,0.5), "lines")
                )
g2 <- g2 + theme(legend.direction="horizontal",
                 legend.position="bottom",
                 legend.justification=c(0,0),
                 legend.box="vertical",
                 legend.box.just="left",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 plot.margin=unit(c(-2,1,0.5,0.5), "lines")
                )
g3 <- arrangeGrob(g1, top = title.A)
g4 <- arrangeGrob(g2, top = title.B)

pdf("NEW_POOLED_Alpha.pdf", height=8,width=8)
lay=rbind(c(1,1),c(2,2))
grid.arrange(g3,g4,layout_matrix=lay)
dev.off()

##### BETA Diversity ######
library(GUniFrac)
library(ape)

phylist <- lapply(seq(1,2),function(i) phyloseq(otu_table(counts(Ldds[[i]],normalize=T),taxa_are_rows=T),sample_data(as.data.frame(colData[[i]])),tax_table(tax_table(myfiltbiom[[i]]))))

phylist <- lapply(phylist,function(o) prune_taxa(rowSums(otu_table(o))>5,o))

nj.ITS
nj.16S

mytree <- list(Bacteria=phy_tree(nj.16S),Fungi=phy_tree(nj.ITS))
mytree <- lapply(mytree,function(o) root(o,outgroup=200,resolve.root=T))

unifracs <- lapply(seq(1,2),function(i) GUniFrac(t(as.matrix(otu_table(phylist[[i]]))),mytree[[i]],alpha=c(0, 0.5, 1))$unifracs)

dw_16 <- unifracs[, , "d_1"] # Weighted UniFrac
du_16 <- unifracs[, , "d_UW"] # Unweighted UniFrac

mytree <- phy_tree(nj.ITS)
mytree <- root(mytree,outgroup=200,resolve.root=T)

unifracs <- GUniFrac(t(as.matrix(otu_table(phylist[[2]]))),mytree,alpha=c(0, 0.5, 1))$unifracs

dw <- lapply(unifracs,function(o) o[, , "d_1"]) # Weighted UniFrac
du <- lapply(unifracs,function(o) o[, , "d_UW"]) # Unweighted UniFrac

lapply(seq(1,2),function(i) adonis(as.dist(du[[i]])~location+Orchard*Sample,colData[[i]],parallel=12,permutations=9999))
lapply(seq(1,2),function(i) adonis(as.dist(dw[[i]])~location+Orchard*Sample,colData[[i]],parallel=12,permutations=9999))

heatPlot <- function(m,s="A") {
  g <- plotHeatmap(m,textSize=15)
  g <- g +theme(axis.title=element_blank(),plot.margin = unit(c(1,1,3,3), "lines"))
  g <- g + annotation_custom(textGrob("Grass",gp = gpar(fontsize = 15)),xmin=11.92,xmax=11.92,ymin=-4,ymax=-4)
  g <- g + annotation_custom(textGrob("Tree",gp = gpar(fontsize = 15)),xmin=35.75,xmax=35.75,ymin=-4,ymax=-4)
  g <- g + annotation_custom(textGrob("Grass",gp = gpar(fontsize = 15)),xmin=59.58,xmax=59.58,ymin=-4,ymax=-4)
  g <- g + annotation_custom(textGrob("Tree",gp = gpar(fontsize = 15)),xmin=83.42,xmax=83.42,ymin=-4,ymax=-4)
  g <- g + annotation_custom(textGrob("\nCider",gp = gpar(fontsize = 15)),xmin=23.83,xmax=23.83,ymin=-6,ymax=-6)
  g <- g + annotation_custom(textGrob("\nDessert",gp = gpar(fontsize = 15)),xmin=71.42,xmax=71.42,ymin=-6,ymax=-6)

  g <- g + annotation_custom(textGrob("\nGrass",rot=90,gp = gpar(fontsize = 15)),ymin=11.92,ymax=11.92,xmin=-9,xmax=-9)
  g <- g + annotation_custom(textGrob("\nTree",rot=90,gp = gpar(fontsize = 15)),ymin=35.75,ymax=35.75,xmin=-9,xmax=-9)
  g <- g + annotation_custom(textGrob("\nGrass",rot=90,gp = gpar(fontsize = 15)),ymin=59.58,ymax=59.58,xmin=-9,xmax=-9)
  g <- g + annotation_custom(textGrob("\nTree",rot=90,gp = gpar(fontsize = 15)),ymin=83.42,ymax=83.42,xmin=-9,xmax=-9)
  g <- g + annotation_custom(textGrob("Cider",rot=90,gp = gpar(fontsize = 15)),ymin=23.83,ymax=23.83,xmin=-11,xmax=-11)
  g <- g + annotation_custom(textGrob("Dessert",rot=90,gp = gpar(fontsize = 15)),ymin=71.42,ymax=71.42,xmin=-11,xmax=-11)

  g <- g + annotation_custom(linesGrob(), xmin = 23.83, xmax = 23.83, ymin = 0, ymax = -1)
  g <- g + annotation_custom(linesGrob(), xmin = 47.75, xmax = 47.75, ymin = 0, ymax = -4)
  g <- g + annotation_custom(linesGrob(), xmin = 71.42, xmax = 71.42, ymin = 0, ymax = -1)

  g <- g + annotation_custom(linesGrob(), ymin = 23.83, ymax = 23.83, xmin = 0, xmax = -1)
  g <- g + annotation_custom(linesGrob(), ymin = 47.75, ymax = 47.75, xmin = 0, xmax = -4)
  g <- g + annotation_custom(linesGrob(), ymin = 71.42, ymax = 71.42, xmin = 0, xmax = -1) 

  g <- g + annotation_custom(textGrob(s,gp = gpar(fontsize = 19)),xmin=100,xmax=100,ymin=95,ymax=95)
  g <- ggplot_gtable(ggplot_build(g))
  g$layout$clip[g$layout$name == "panel"] <- "off"
  return(g)
}

pdf("POOLED_betaX.pdf",height=9,width=11)
grid.arrange(
  heatPlot(du[[2]],"A"),
  heatPlot(dw[[2]],"B"),
  heatPlot(du[[1]],"C"),
  heatPlot(dw[[1]],"D"),
  nrow=2,
  ncol=2
)
dev.off()




#===============================================================================
#       differential analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(	
	(sample_data(obj)[[10]]!="duplicate")&
	(sample_data(obj)[[1]]!="C")&
	(colSums(otu_table(obj))>999)
,obj))


myfiltbioms <- lapply(myfiltbioms,function(obj) {
	obj@sam_data$location <- as.factor(obj@sam_data$meters)
	colnames(sample_data(obj))[c(1,6,11)] <- c("Sample","Distance","Orchard")
	levels(sample_data(obj)[[1]]) <- c("Aisle","Tree")
	return(obj)
})


Ldds <- lapply(myfiltbioms,function(obj) phylo_to_des(obj,fit=F,calcFactors=geoSet))

Ldds$Fungi$Group <- as.factor(paste(Ldds$Fungi$Orchard,Ldds$Fungi$Sample,Ldds$Fungi$location,sep="_"))
Ldds$Bacteria$Group <- as.factor(paste(Ldds$Bacteria$Orchard,Ldds$Bacteria$Sample,Ldds$Bacteria$location,sep="_"))
Ldds<-lapply(Ldds,function(o) collapseReplicates2(o,groupby=o$Group))
Ldds$Fungi$Orchard<- droplevels(Ldds$Fungi$Orchard)
Ldds$Fungi$Sample<- droplevels(Ldds$Fungi$Sample)
Ldds$Fungi$location<- droplevels(Ldds$Fungi$location)
Ldds$Bacteria$Orchard<- droplevels(Ldds$Bacteria$Orchard)
Ldds$Bacteria$Sample<- droplevels(Ldds$Bacteria$Sample)
Ldds$Bacteria$location<- droplevels(Ldds$Bacteria$location)


#Ldds <- lapply(Ldds,function(obj) obj[rowSums(counts(obj)>5)>=3,])

design=~Orchard + Sample + Orchard:location + Orchard:Sample

lapply(seq(1,length(Ldds)),function(i) design(Ldds[[i]])<<-design)

Ldds <- lapply(Ldds,function(obj) DESeq(obj,parallel=T))

alpha <- 0.05

contrast=c("Sample","Aisle", "Tree" ) # main effect
contrast=c("Orchard","Cider","Dessert") # orchard effect
contrast=list( "OrchardCider.SampleAisle","OrchardCider.SampleTree") # Grass:Tree (cider)
contrast=list( "OrchardDessert.SampleAisle","OrchardDessert.SampleTree") # Grass:Tree (Dessert)

res <- lapply(Ldds,function(obj) results(obj,alpha=alpha,parallel=T,contrast=contrast))

#res.merge <- lapply(res,function(o) lapply(seq(1:2),function(i) 

res.merge <- lapply(seq(1,length(Ldds)),function(i) 
	data.table(inner_join(
		data.table(OTU=rownames(res[[i]]),as.data.frame(res[[i]])),
		data.table(OTU=rownames(tax_table(mybioms[[i]])),as.data.frame(as.matrix(tax_table(mybioms[[i]])))) 
	))
)


write.table(res.merge[[1]],"Pooled_bacteria.dessert.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[2]],"Pooled_fungi.dessert.res",sep="\t",quote=F,na="",row.names=F)



## figure S5 

t1 <- lapply(res.main,function(o) countTaxa2(taxaConf(o[o$padj<=0.05&o$log2FoldChange<0,8:21],0.65,3),"rank"))
lapply(seq(1,2),function(i) colnames(t1[[i]])[ncol(t1[[i]])]<<-"Tree")
t1<- lapply(seq(1,2),function(i) full_join(t1[[i]],countTaxa2(taxaConf(res.main[[i]][res.main[[i]]$padj<=0.05&res.main[[i]]$log2FoldChange>0,8:21],0.65,3),"rank"),by="Taxa"))
lapply(seq(1,2),function(i) colnames(t1[[i]])[ncol(t1[[i]])]<<-"Grass")
t1 <- lapply(seq(1,2),function(i) full_join(t1[[i]],countTaxa2(taxaConf(res.orchard[[i]][res.orchard[[i]]$padj<=0.05,8:21],0.65,3),"rank"),by="Taxa"))
lapply(seq(1,2),function(i) colnames(t1[[i]])[ncol(t1[[i]])]<<-"Orchard")
t1 <- lapply(seq(1,2),function(i) full_join(t1[[i]],countTaxa2(taxaConf(res.cider[[i]][res.cider[[i]]$padj<=0.05,8:21],0.65,3),"rank"),by="Taxa"))
lapply(seq(1,2),function(i) colnames(t1[[i]])[ncol(t1[[i]])]<<-"Interaction")

names(t1) <- c("Bacteria","Fungi")
lapply(seq(1,2),function(i) write.table(t1[[i]],paste("POOLED_",names(t1)[i],".sig.class.count.txt",sep=""), sep="\t", quote=F,na="0",row.names=F))

ggcorr(t1[[1]][,2:5],label=T,label_round=2,label_size=6,legend.size=12,size=6)
ggcorr(t1[[2]][,2:5],label=T,label_round=2,label_size=6,legend.size=12,size=6)


## MA plots

pdf("POOLED_ma.plots.pdf")
lapply(res.main,function(o) plot_ma(o,xlims=c(-10,10),textsize=20))
lapply(res.orchard,function(o) plot_ma(o,xlims=c(-10,10),textsize=20))
lapply(res.cider,function(o) plot_ma(o,xlims=c(-10,10),textsize=20))
lapply(res.dessert,function(o) plot_ma(o,xlims=c(-10,10),textsize=20))
dev.off()




 ggpairs(test,aes(colour=colour,alpha=0.4),
	upper = list(continuous = "cor", combo = "facetdensity"),
	lower = list(continuous = "smooth", combo = "dot_no_facet"),
	diag = list(continuous="densityDiag"))
dev.off()

g<- ggpairs(test,aes(colour=colour,alpha=0.4))
g <- g + scale_colour_manual(values=cbbPalette)
g <- g + theme_bw()
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(size=16))


	obj<-res.merge[[1]]
	with(obj,plot(y=log2FoldChange,x=log10(baseMean),pch=20, xlim=c(-6,6),bty="n", cex.lab=1.5, cex.axis=1.4, cex.main=1.4, cex.sub=1.4,
	xlab=expression("Log"[2]*" Fold Change"),ylab=expression("Log"[10]*" Mean Expression")))
	with(subset(obj, padj<0.05 ), points(y=log2FoldChange, x=log10(baseMean), pch=20, col="#E69F00"))
	with(subset(obj, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="#56B4E9"))
	with(subset(obj, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="#009E73"))


