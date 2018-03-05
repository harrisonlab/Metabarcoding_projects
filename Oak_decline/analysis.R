#===============================================================================
#       Load libraries
#===============================================================================
library(DESeq2)
library(BiocParallel)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
register(MulticoreParam(12))
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)


#===============================================================================
#       Load data 
#===============================================================================

##### 16S #####

#biom_file = "16S.taxa.biom"
#colData = "colData"
#mybiom <- import_biom(biom_file) 
#sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
#tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65)
#biom16<-mybiom

##### ITS #####

#biom_file = "ITS.taxa.biom"
#colData = "colData"
#mybiom <- import_biom(biom_file) 
#sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
#tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65)
#biomITS<-mybiom
#mybioms <- list(bacteria=biom16,fungi=biomITS)

ubiom_BAC <- loadData("16S.otu_table.txt","colData","16S.taxa","16S.phy",RHB="BAC")
ubiom_FUN <- loadData("ITS.otu_table.txt","colData","ITS.taxa","ITS.phy",RHB="FUN")

#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("ITS.taxa")
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData

#===============================================================================
#       Filter data
#===============================================================================

colData <- colData[colData$OAK==1,]
dds <- dds[rowSums(counts(dds, normalize=T))>4,]


#===============================================================================
#       PCA analysis
#===============================================================================

##### pca all sites #####

myfiltbioms <- lapply(mybioms,function(obj) prune_taxa(rowSums(otu_table(obj))>5,obj))
mypcas <-  lapply(myfiltbioms, function(obj) plotPCA(obj,design="1",returnData=T,calcFactors=geoSet))

dfs <-lapply(seq(1:2), function(i) t(data.frame(t(mypcas[[i]]$x)*mypcas[[i]]$percentVar)))

pdf("all_pca2.pdf")
lapply(seq(1:2),function(i) plotOrd(dfs[[i]],sample_data(myfiltbioms[[i]]),design=c("site","condition"),dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)+ ggtitle(names(myfiltbioms)[i]))
lapply(seq(1:2),function(i) plotOrd(dfs[[i]],sample_data(myfiltbioms[[i]]),design=c("site","condition"),dimx=2,dimy=3,xlabel="PC2",ylabel="PC3",pointSize=1.5,cbPallete=T)+ ggtitle(names(myfiltbioms)[i]))
dev.off()

###### pca chestnuts/bigwood #####

myfiltbioms <- lapply(mybioms,function(obj) prune_taxa(rowSums(otu_table(obj))>5,obj))

filterFun=function(o,f){
	o<-prune_samples(sample_data(o)$site==f,o)
}

mypcasites <- lapply(myfiltbioms,function(obj) 
  sapply(seq(1:2),function(i) 
    setNames(
      list(
        plotPCA(obj,design="1",returnData=T,calcFactors=geoSet,filterFun=filterFun,filter=levels(sample_data(obj)$site)[i])
      ),c(levels(sample_data(obj)$site)[i])
    )
  )
)

dfs <- lapply(mypcasites, function(obj)
	lapply(obj, function(o) t(data.frame(t(o$x)*o$percentVar)))
)

pdf("site_pca_all_2.pdf")
glist1 <- lapply(dfs,function(obj) lapply(obj,function(o)
	plotOrd(o,sample_data(biom16)[rownames(sample_data(biom16))%in%rownames(o)],design="condition",dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
))
glist2 <- lapply(dfs,function(obj) lapply(obj,function(o)
	plotOrd(o,sample_data(biom16)[rownames(sample_data(biom16))%in%rownames(o)],design="condition",dimx=2,dimy=3,xlabel="PC2",ylabel="PC3",pointSize=1.5,cbPallete=T)
))
nn <-  sapply(names(glist1),function(x) paste(x,sapply(glist1,names)))
lapply(rep(1:2),function(i) lapply(rep(1:2),function(j) glist1[[i]][[j]]+ggtitle(nn[j,i])))
lapply(rep(1:2),function(i) lapply(rep(1:2),function(j) glist2[[i]][[j]]+ggtitle(nn[j,i])))
dev.off()


##### pca chestnuts - with location info #####

mypca_chestnuts <- list(bacteria=mypcasites[[1]][[2]],fungi=mypcasites[[2]][[2]])


#===============================================================================
#       differential analysis
#===============================================================================

Ldds <- lapply(mybioms,function(o) 
  sapply(seq(1:2),function(i) 
    setNames(
      list(phylo_to_des(filterFun(o,levels(sample_data(o)$site)[i]),fit=F,calcFactors=geoSet)),
      c(levels(sample_data(o)$site)[i])
    )
  )
)

Ldds <- lapply(Ldds,function(obj) lapply(obj,function(o) o[rowSums(counts(o))>5,]))
designs=c(formula(~condition),formula(~tree+condition))
lapply(seq(1:2),function(i) lapply(seq(1:2), function(j) design(Ldds[[i]][[j]]) <<- designs[[j]])) # note the <<- for global assignment 
Ldds <- lapply(Ldds,function(obj) lapply(obj,function(o) DESeq(o,parallel=T)))

alpha <- 0.05

res <- lapply(Ldds,function(obj)lapply(obj,function(o)results(o,alpha=alpha,parallel=T)))

#res.merge <- lapply(res,function(o) lapply(seq(1:2),function(i) 

res.merge <- lapply(seq(1:2),function(i) lapply(seq(1:2),function(j) 
	data.table(inner_join(
		data.table(OTU=rownames(res[[i]][[j]]),as.data.frame(res[[i]][[j]])),
		data.table(OTU=rownames(tax_table(mybioms[[i]])),as.data.frame(as.matrix(tax_table(mybioms[[i]]))))
	))
))

write.table(res.merge[[1]][[1]],"Bigwood.bacteria.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[2]][[1]],"Bigwood.fungi.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[1]][[2]],"Chestnuts.bacteria.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[2]][[2]],"Chestnuts.fungi.res",sep="\t",quote=F,na="",row.names=F)

## Volcano plots

pdf("dispersion_plots.pdf")

lapply(unlist(res.merge,recursive=F),function(o) {
	with(o,plot(log2FoldChange,log10(baseMean),pch=20, xlim=c(-2.5,2)))
	with(subset(o, padj<.05 ), points(log2FoldChange, log10(baseMean), pch=20, col="red"))
	with(subset(o, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="orange"))
	with(subset(o, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="green"))
})

dev.off()



#===============================================================================
#       Alpha diversity analysis
#===============================================================================

pdf("16S_v2.alpha.pdf")
plot_richness(biom16,color="condition",x="site",measures=c("Chao1", "Shannon", "Simpson"))
dev.off()
pdf("ITS_v2.alpha.pdf")
plot_richness(biomITS,color="condition",x="site",measures=c("Chao1", "Shannon", "Simpson"))
dev.off()

all_alpha <- plot_richness(biom16,returnData=T)


summary(aov(Chao1~condition*site,all_alpha))[[1]][2]
summary(aov(Shannon~location+(Sample*Orchard),all_alpha))[[1]][2]
summary(aov(Simpson~location+(Sample*Orchard),all_alpha))[[1]][2]


#===============================================================================
#       Beta diversity analysis
#===============================================================================

library(ape)
library(data.table)
library(future)
plan(multiprocess)

phylip_data_ITS = fread("ITS.phy",skip=1)
ftITS <- future({nj(as.dist(data.frame(phylip_data_ITS,row.names="V1")))})
nj.ITS <- value(ftITS) 
write.tree(nj.ITS,"ITS.tree")
phy_tree(biomITS) <- nj.ITS


phylip_data_16S = fread("16S.phy",skip=1)
temp <- as.dist(data.frame(phylip_data_16S,row.names="V1"))
ft <- future({nj(temp)})
nj.16S <- value(ft) 
write.tree(nj.16S,"ITS.tree")
phy_tree(biom16S) <- nj.16S



#===============================================================================
#       network analysis
#===============================================================================

library(cooccur)
library(future) 
plan(multiprocess) 


myfiltbiom <- prune_samples(sample_data(biom16)$site=="Chestnuts",biom16)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 


cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1 


fht16 <- future({cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)})
fst16 <- future({cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)})
#fs16 <- future({cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T))}

CHcoHmodel16 <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)
CHcoSmodel16 <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)
 
 
CHcoHmodel16$results$padj <- p.adjust(apply(CHcoHmodel16$results[,8:9],1, min),"BH") 
CHcoSmodel16$results$padj <- p.adjust(apply(CHcoSmodel16$results[,8:9],1, min),"BH") 
#CHcomodel$16results$padj <- p.adjust(apply(CHcomodel$results[,8:9],1, min),"BH") 


myfiltbiom <- prune_samples(sample_data(biom16)$site=="Bigwood",biom16) 
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom))) 
 
cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Control",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 

cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1  

BWfht16 <- future({cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)})
BWfst16 <- future({cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)})

BWcoHmodel16 <- value(BWfht16) #  cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T) #value(BWfth1) 
BWcoSmodel16 <- value(BWfst16) #  cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T) #value(BWfts1) 
#BWcomodel16 <-  cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T) 

BWcoHmodel16$results$padj <- p.adjust(apply(BWcoHmodel16$results[,8:9],1, min),"BH") 
BWcoSmodel16$results$padj <- p.adjust(apply(BWcoSmodel16$results[,8:9],1, min),"BH") 

##Speculation

myfiltbiom <- prune_samples(sample_data(biomITS)$site=="Speculation",biomITS)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))
cotable <- cotable[rowSums(cotable)>5,colSums(cotable)>5]

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 

cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1  


SPcoHmodel <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)
SPcoSmodel <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)


SPcoHmodel $results$padj <- p.adjust(apply(SPcoHmodel $results[,8:9],1, min),"BH") 
SPcoSmodel $results$padj <- p.adjust(apply(SPcoSmodel $results[,8:9],1, min),"BH") 


cotable[cotable>0] <- 1 

SPcomodel <- cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T)
SPcomodel$results$padj <- p.adjust(apply(SPcomodel$results[,8:9],1, min),"BH") 


myfiltbiom <- prune_samples(sample_data(biom16)$site=="Speculation",biom16)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 

cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1  


SPcoHmodel16 <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)
SPcoSmodel16 <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)

SPcoHmodel16 $results$padj <- p.adjust(apply(SPcoHmodel16 $results[,8:9],1, min),"BH") 
SPcoSmodel16 $results$padj <- p.adjust(apply(SPcoSmodel16 $results[,8:9],1, min),"BH") 

cotable <- cotable[rowSums(cotable)>5,colSums(cotable)>5]
cotable[cotable>0] <- 1 

SPcomodel16 <- cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T)
SPcomodel16$results$padj <- p.adjust(apply(SPcomodel16$results[,8:9],1, min),"BH") 

## ITS
summary(CHcoHmodel)
summary(CHcoSmodel) 
summary(BWcoHmodel)
summary(BWcoSmodel)
summary(SPcomodel)

nrow(CHcoHmodel$results[CHcoHmodel$results$padj<=0.1&CHcoHmodel$results$p_gt<=0.1,])
nrow(CHcoHmodel$results[CHcoHmodel$results$padj<=0.1&CHcoHmodel$results$p_lt<=0.1,])
nrow(CHcoSmodel$results[CHcoSmodel$results$padj<=0.1&CHcoSmodel$results$p_gt<=0.1,])
nrow(CHcoSmodel$results[CHcoSmodel$results$padj<=0.1&CHcoSmodel$results$p_lt<=0.1,])
nrow(BWcoHmodel$results[BWcoHmodel$results$padj<=0.1&BWcoHmodel$results$p_gt<=0.1,])
nrow(BWcoHmodel$results[BWcoHmodel$results$padj<=0.1&BWcoHmodel$results$p_lt<=0.1,])
nrow(BWcoSmodel$results[BWcoSmodel$results$padj<=0.1&BWcoSmodel$results$p_gt<=0.1,])
nrow(BWcoSmodel$results[BWcoSmodel$results$padj<=0.1&BWcoSmodel$results$p_lt<=0.1,])
nrow(SPcoHmodel$results[SPcoHmodel$results$padj<=0.1&SPcoHmodel$results$p_gt<=0.1,])
nrow(SPcoHmodel$results[SPcoHmodel$results$padj<=0.1&SPcoHmodel$results$p_lt<=0.1,])
nrow(SPcoSmodel$results[SPcoSmodel$results$padj<=0.1&SPcoSmodel$results$p_gt<=0.1,])
nrow(SPcoSmodel$results[SPcoSmodel$results$padj<=0.1&SPcoSmodel$results$p_lt<=0.1,])


pdf("obs.exp.fungi.pdf")
obs.v.exp(CHcoHmodel)
obs.v.exp(CHcoHmodel)
obs.v.exp(BWcoHmodel)
obs.v.exp(BWcoSmodel)
obs.v.exp(SPcoHmodel)
obs.v.exp(SPcoSmodel)
dev.off()


## 16S

summary(CHcoHmodel16)
summary(CHcoSmodel16)
summary(BWcoHmodel16)
summary(BWcoSmodel16)
summary(SPcomodel16)

nrow(CHcoHmodel16$results[CHcoHmodel16$results$padj<=0.1&CHcoHmodel16$results$p_gt<=0.1,])
nrow(CHcoHmodel16$results[CHcoHmodel16$results$padj<=0.1&CHcoHmodel16$results$p_lt<=0.1,])
nrow(CHcoSmodel16$results[CHcoSmodel16$results$padj<=0.1&CHcoSmodel16$results$p_gt<=0.1,])
nrow(CHcoSmodel16$results[CHcoSmodel16$results$padj<=0.1&CHcoSmodel16$results$p_lt<=0.1,])
nrow(BWcoHmodel16$results[BWcoHmodel16$results$padj<=0.1&BWcoHmodel16$results$p_gt<=0.1,])
nrow(BWcoHmodel16$results[BWcoHmodel16$results$padj<=0.1&BWcoHmodel16$results$p_lt<=0.1,])
nrow(BWcoSmodel16$results[BWcoSmodel16$results$padj<=0.1&BWcoSmodel16$results$p_gt<=0.1,])
nrow(BWcoSmodel16$results[BWcoSmodel16$results$padj<=0.1&BWcoSmodel16$results$p_lt<=0.1,])
nrow(SPcoHmodel16$results[SPcoHmodel16$results$padj<=0.1&SPcoHmodel16$results$p_gt<=0.1,])
nrow(SPcoHmodel16$results[SPcoHmodel16$results$padj<=0.1&SPcoHmodel16$results$p_lt<=0.1,])
nrow(SPcoSmodel16$results[SPcoSmodel16$results$padj<=0.1&SPcoSmodel16$results$p_gt<=0.1,])
nrow(SPcoSmodel16$results[SPcoSmodel16$results$padj<=0.1&SPcoSmodel16$results$p_lt<=0.1,])


head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
head(CHcoHmodel$results[order(CHcoHmodel$results$p_gt,CHcoHmodel$results$padj),])
head(CHcoSmodel$results[order(CHcoSmodel$results$p_lt,CHcoSmodel$results$padj),])
head(CHcoSmodel$results[order(CHcoSmodel$results$p_gt,CHcoSmodel$results$padj),])
head(SPcoHmodel$results[order(SPcoHmodel$results$p_lt,SPcoHmodel$results$padj),])
head(SPcoHmodel$results[order(SPcoHmodel$results$p_gt,SPcoHmodel$results$padj),])
head(SPcoSmodel$results[order(SPcoSmodel$results$p_lt,SPcoSmodel$results$padj),])
head(SPcoSmodel$results[order(SPcoSmodel$results$p_gt,SPcoSmodel$results$padj),])
head(BWcoHmodel$results[order(BWcoHmodel$results$p_lt,BWcoHmodel$results$padj),])
head(BWcoHmodel$results[order(BWcoHmodel$results$p_gt,BWcoHmodel$results$padj),])
head(BWcoSmodel$results[order(BWcoSmodel$results$p_lt,BWcoSmodel$results$padj),])
head(BWcoSmodel$results[order(BWcoSmodel$results$p_gt,BWcoSmodel$results$padj),])

X <- head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_lt,CHcoHmodel16$results$padj),])
X <- rbind(X,head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_gt,CHcoHmodel16$results$padj),]))
X <- rbind(X,head(CHcoSmodel16$results[order(CHcoSmodel16$results$p_lt,CHcoSmodel16$results$padj),]))
X <- rbind(X,head(CHcoSmodel16$results[order(CHcoSmodel16$results$p_gt,CHcoSmodel16$results$padj),]))
X <- rbind(X,head(SPcoHmodel16$results[order(SPcoHmodel16$results$p_lt,SPcoHmodel16$results$padj),]))
X <- rbind(X,head(SPcoHmodel16$results[order(SPcoHmodel16$results$p_gt,SPcoHmodel16$results$padj),]))
X <- rbind(X,head(SPcoSmodel16$results[order(SPcoSmodel16$results$p_lt,SPcoSmodel16$results$padj),]))
X <- rbind(X,head(SPcoSmodel16$results[order(SPcoSmodel16$results$p_gt,SPcoSmodel16$results$padj),]))
X <- rbind(X,head(BWcoHmodel16$results[order(BWcoHmodel16$results$p_lt,BWcoHmodel16$results$padj),]))
X <- rbind(X,head(BWcoHmodel16$results[order(BWcoHmodel16$results$p_gt,BWcoHmodel16$results$padj),]))
X <- rbind(X,head(BWcoSmodel16$results[order(BWcoSmodel16$results$p_lt,BWcoSmodel16$results$padj),]))
X <- rbind(X,head(BWcoSmodel16$results[order(BWcoSmodel16$results$p_gt,BWcoSmodel16$results$padj),]))


Y <- head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
Y <- rbind(Y,head(CHcoHmodel$results[order(CHcoHmodel$results$p_gt,CHcoHmodel$results$padj),]))
Y <- rbind(Y,head(CHcoSmodel$results[order(CHcoSmodel$results$p_lt,CHcoSmodel$results$padj),]))
Y <- rbind(Y,head(CHcoSmodel$results[order(CHcoSmodel$results$p_gt,CHcoSmodel$results$padj),]))
Y <- rbind(Y,head(SPcoHmodel$results[order(SPcoHmodel$results$p_lt,SPcoHmodel$results$padj),]))
Y <- rbind(Y,head(SPcoHmodel$results[order(SPcoHmodel$results$p_gt,SPcoHmodel$results$padj),]))
Y <- rbind(Y,head(SPcoSmodel$results[order(SPcoSmodel$results$p_lt,SPcoSmodel$results$padj),]))
Y <- rbind(Y,head(SPcoSmodel$results[order(SPcoSmodel$results$p_gt,SPcoSmodel$results$padj),]))
Y <- rbind(Y,head(BWcoHmodel$results[order(BWcoHmodel$results$p_lt,BWcoHmodel$results$padj),]))
Y <- rbind(Y,head(BWcoHmodel$results[order(BWcoHmodel$results$p_gt,BWcoHmodel$results$padj),]))
Y <- rbind(Y,head(BWcoSmodel$results[order(BWcoSmodel$results$p_lt,BWcoSmodel$results$padj),]))
Y <- rbind(Y,head(BWcoSmodel$results[order(BWcoSmodel$results$p_gt,BWcoSmodel$results$padj),]))



#===============================================================================
#       Speculation site
#===============================================================================

# PCA

filterFun=function(o,f){
	o<-prune_samples(sample_data(o)[[2]]=="Speculation",o)
	prune_taxa(rowSums(otu_table(o))>5,o)
}

mypca <- plotPCA(biomITS,design="1",returnData=T,filterFun=filterFun,calcFactors=geoMeans)
myfiltbiom <- prune_samples(sample_data(biomITS)[[2]]=="Speculation",biomITS)
myfiltbiom <-  prune_taxa(rowSums(otu_table(myfiltbiom ))>5,myfiltbiom )
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

pdf("spec_its.pdf")
plotOrd(df,sample_data(myfiltbiom),design="condition",dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
dev.off

mypca <- plotPCA(biom,design="1",returnData=T,filterFun=filterFun,calcFactors=geoMeans)
myfiltbiom <- prune_samples(sample_data(biom)[[2]]=="Speculation",biom)
myfiltbiom <-  prune_taxa(rowSums(otu_table(myfiltbiom ))>5,myfiltbiom )
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

pdf("spec_its.pdf")
plotOrd(df,sample_data(myfiltbiom),design="condition",dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
dev.off()


# differential analysis

site="Speculation"
dds <- phylo_to_des(biomITS ,fit=F,calcFactors=geoMeans)
dds <- dds[ rowSums(counts(dds)) > 5, ]
#dds <- dds[,dds$site==site]
dds$condition <- droplevels(dds$condition)
dds$site <- droplevels(dds$site)
dds$tree <- droplevels(dds$tree)

design=~condition#~tree + condition
design(dds) <- design
dds <- DESeq(dds,parallel=T)
alpha <- 0.05

contrast=c("condition","Healthy","Symptom")
res <- results(dds,contrast=contrast,alpha=alpha,parallel=T)

myfiltbiom <- prune_samples(sample_data(biomITS)$site==site,mybiomITS)

res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)

write.table(res.merge,"speculation.fun.res",sep="\t",quote=F,na="",row.names=F)


dds <- phylo_to_des(biom ,fit=F)
dds <- dds[ rowSums(counts(dds)) > 5, ]
#dds <- dds[,dds$site==site]
dds$condition <- droplevels(dds$condition)
dds$site <- droplevels(dds$site)
dds$tree <- droplevels(dds$tree)

design=~condition#~tree + condition
design(dds) <- design
dds <- DESeq(dds,parallel=T)
alpha <- 0.05

contrast=c("condition","Healthy","Symptom")
res <- results(dds,contrast=contrast,alpha=alpha,parallel=T)

myfiltbiom <- prune_samples(sample_data(biom)$site==site,biom)

res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)

write.table(res.merge,"speculation.bac.res",sep="\t",quote=F,na="",row.names=T)


## Volcano plots

pdf("volcano_spec.pdf")

with(res.merge.its,plot(log2FoldChange,log10(baseMean),pch=20, xlim=c(-2.5,2)))
with(subset(res.merge.its, padj<0.05 ), points(log2FoldChange, log10(baseMean), pch=20, col="red"))
with(subset(res.merge.its, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="orange"))
with(subset(res.merge.its, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="green"))


dev.off()






### Other stuff


# using combined biological replicates (I think this is a bad idea..)
groupby=paste(dds$condition,dds$tree,sep="_")
dds2 <- collapseReplicates2(dds,groupby=groupby) # uses mean for grouping rather than sum
design(dds2) <- ~condition
dds2 <- DESeq(dds2,parallel=T)
res2 <- results(dds2,alpha=alpha,parallel=T)

pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$tree))

df <- t(data.frame(t(mypca$x)*mypca$percentVar))
d <- t(data.frame(t(pc.res)*mypca$percentVar))

sample_data(myfiltbiom)$loc <- as.numeric(levels(sample_data(myfiltbiom)$tree)[sample_data(myfiltbiom)$tree])
plotOrd(df,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=1,dimy=2,xlabel="PC1",ylabel="PC2")
plotOrd(d,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=1,dimy=2,xlabel="PC1",ylabel="PC2")
plotOrd(df,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=2,dimy=3,xlabel="PC2",ylabel="PC3")
plotOrd(d,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=2,dimy=3,xlabel="PC2",ylabel="PC3")


pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$tree))

mybiom <- biom
#mybiom <- biomITS

myfiltbiom <-  prune_taxa(rowSums(otu_table(mybiom ))>5,mybiom )
mypca <- plotPCA(myfiltbiom ,design="1",returnData=T,calcFactors=geoMeans,filterFun=filterFun,filter=...)


#mypca <- plotPCA(myfiltbiom ,design="1",returnData=T,calcFactors=geoMeans) # ITS only (bypass log geometric means error)


df <- t(data.frame(t(mypca$x)*mypca$percentVar))
pdf("all_S_pca.pdf")
#pdf("all_ITS_pca.pdf")

plotOrd(df,sample_data(myfiltbiom),design=c("site","condition"),dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
plotOrd(df,sample_data(myfiltbiom),design=c("site","condition"),dimx=2,dimy=3,xlabel="PC2",ylabel="PC3",pointSize=1.5,cbPallete=T)
dev.off()
