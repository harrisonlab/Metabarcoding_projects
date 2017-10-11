#===============================================================================
#       Load libraries
#===============================================================================

# library(phyloseq)
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
# library(cooccur)
# library(parallel)

#===============================================================================
#       Load data 
#===============================================================================

# load denoised otu count table
countData <- read.table("BAC.zotus_table.txt",header=T,sep="\t",row.names=1, comment.char = "")

# load sample metadata
colData <- read.table("colData",header=T,sep="\t",row.names=1)

# load taxonomy data
taxData <- read.table("zBAC.taxa",header=F,sep=",",row.names=1)

# reorder columns
taxData<-taxData[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]

# add best "rank" at 0.65 confidence and tidy-up the table
taxData<-phyloTaxaTidy(taxData,0.65)

# save data into a list
ubiom_BAC <- list(
	countData=countData,
	colData=colData,
	taxData=taxData,
	RHB="BAC"
)

# Fungi all in one call
ubiom_FUN <- list(
	countData=read.table("FUN.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=T,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 

# Oomycetes
ubiom_OO <- list(
	countData=read.table("OO.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData2",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zOO.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="OO"
) 

# Nematodes 
ubiom_NEM <- list(
	countData=read.table("NEM.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData2",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zNEM.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="OO"
) 

#===============================================================================
#       Combine species 
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species). Works well for Oomycetes and Fungi

# attach OO objects
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
# list of species with more than one associated OTU
combinedTaxa <- combineTaxa("zOO.taxa")
# show the list
combinedTaxa[,1]
# manual filter list to remove none species (e.g. unknown, Pythium aff)
combinedTaxa <- combinedTaxa[c(-2,-8,-9,-17,-20,-25),]
# adjust countData for combined taxa
countData <- combCounts(combinedTaxa,countData)
# adjust taxData for combined taxa
taxData <- combTaxa(combinedTaxa,taxData)
# copy back to object
ubiom_OO$countData <- countData
ubiom_OO$taxData <- taxData

# list of species with more than one associated OTU
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("zFUN.taxa")
combinedTaxa[,1]
combinedTaxa <- combinedTaxa
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData



#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# remove low count samples and control samples (not needed here)
filter <- (colSums(counts(dds))>=1000) & (colData$pair!="C")
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - use geoMeans function if
# every gene contains at least one zero, as cannot compute log geometric means
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
#sizeFactors(dds) <-geoMeans(dds)
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))


#===============================================================================
#       Filter data 
#============================================================================

# Pythium specific filter to remove OTUs which are unlikely part of the SAR kingdom
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])

dds <- dds[myfilter,]

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot 
colData$location<-as.number(colData$pair)

# plot the PCA
pdf(paste(RHB,"VA.pdf",sep="_"))
plotOrd(df,colData,design="condition",xlabel="PC1",ylabel="PC2")
plotOrd(df,colData,shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()


### remove spatial information (this uses the factor "pair" not the numeric "location") and plot

pc.res <- resid(aov(mypca$x~colData$pair))
ds <- lapply(seq(1,length(pc.res)), function(i) t(data.frame(t(pc.res[[i]])*mypcas[[i]]$percentVar)))


pdf("ND_all_pcas_pairs_noloc.pdf")
lapply(seq(1,length(myfiltbioms)),function(i) plotOrd(ds[[i]],sample_data(myfiltbioms[[i]]),shape="condition",design="location",continuous=T,dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)+ ggtitle(names(myfiltbioms)[i]))
lapply(seq(1,length(myfiltbioms)),function(i) plotOrd(ds[[i]],sample_data(myfiltbioms[[i]]),shape="condition",design="location",continuous=T,dimx=2,dimy=3,xlabel="PC2",ylabel="PC3",pointSize=1.5,cbPallete=T)+ ggtitle(names(myfiltbioms)[i]))
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))
Ldds <- lapply(myfiltbioms,function(obj) phylo_to_des(obj,fit=F,calcFactors=geoSet))

Ldds <- lapply(Ldds,function(obj) obj[rowSums(counts(obj))>5,])
#designs=c(formula(~condition),formula(~pair+condition))

design = ~pair+condition
lapply(seq(1,length(Ldds)),function(i) design(Ldds[[i]])<<-design)

Ldds <- lapply(Ldds,function(obj) DESeq(obj,parallel=T))

alpha <- 0.05


res <- lapply(Ldds,function(obj) results(obj,alpha=alpha,parallel=T,contrast=contrast))

#res.merge <- lapply(res,function(o) lapply(seq(1:2),function(i) 

res.merge <- lapply(seq(1,length(Ldds)),function(i) 
	data.table(inner_join(
		data.table(OTU=rownames(res[[i]]),as.data.frame(res[[i]])),
		data.table(OTU=rownames(tax_table(mybioms[[i]])),as.data.frame(as.matrix(tax_table(mybioms[[i]])))) 
	)))


write.table(res.merge[[1]],"bacteria.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[2]],"fungi.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[3]],"oo.res",sep="\t",quote=F,na="",row.names=F)
write.table(res.merge[[4]],"nem.res",sep="\t",quote=F,na="",row.names=F)

## Volcano plots

pdf("volcano_plots.pdf")
lapply(res.merge,function(obj) {
	with(obj,plot(log2FoldChange,log10(baseMean),pch=20, xlim=c(-5,5)))
	with(subset(obj, padj<.1 ), points(log2FoldChange, log10(baseMean), pch=20, col="red"))
	with(subset(obj, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="orange"))
	with(subset(obj, padj<.1 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="green"))
})
dev.off()


#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()

#===============================================================================
#       network analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

cotables <- lapply(myfiltbioms,function(obj) as.data.frame(as.matrix(otu_table(obj))))
cotables_h <- lapply(seq(1,length(myfiltbioms)),function(i)  cotables[[i]][,row.names(sample_data(prune_samples(sample_data(myfiltbioms[[i]])$condition=="Healthy",myfiltbioms[[i]])))]) 
cotables_s <- lapply(seq(1,length(myfiltbioms)),function(i)  cotables[[i]][,row.names(sample_data(prune_samples(sample_data(myfiltbioms[[i]])$condition=="Symptom",myfiltbioms[[i]])))])

cotables_h <- lapply(cotables_h, function(obj) obj[rowSums(obj)>5,colSums(obj)>5])
cotables_s <- lapply(cotables_s, function(obj) obj[rowSums(obj)>5,colSums(obj)>5])

lapply(seq(1,length(cotables_h)), function(i) cotables_h[[i]][cotables_h[[i]]>0] <<- 1)
lapply(seq(1,length(cotables_s)), function(i) cotables_s[[i]][cotables_s[[i]]>0] <<- 1)

CFcoHmodels <- mclapply(cotables_h, function(obj) cooccur2(obj,type = "spp_site",spp_names = T,thresh = T),mc.cores=4)
CFcoSmodels <- mclapply(cotables_s, function(obj) cooccur2(obj,type = "spp_site",spp_names = T,thresh = T),mc.cores=4)

lapply(seq(1,length(CFcoHmodels)), function(i) {
	CFcoHmodels[[i]]$results$padj <<- p.adjust(apply(CFcoHmodels[[i]]$results[,8:9],1, min),"BH")
	CFcoSmodels[[i]]$results$padj <<- p.adjust(apply(CFcoSmodels[[i]]$results[,8:9],1, min),"BH")
})

lapply(CFcoHmodels, function(obj) {nrow(obj$results[obj$results$padj<=0.1,]})

lapply(CFcoHmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_gt<=0.05,]))
lapply(CFcoHmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_lt<=0.05,]))

lapply(CFcoSmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_gt<=0.05,]))
lapply(CFcoSmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_lt<=0.05,]))


X <- rbind.fill(lapply(CFcoHmodels, function(obj) head(obj$results[order(obj$results$padj),],6)))
Y <- rbind.fill(lapply(CFcoSmodels, function(obj) head(obj$results[order(obj$results$padj),],6)))


HcoHmodel16$results$p_lt

write.table(X,"cooc.healthy.txt",sep="\t",quote=F,row.names=F)

head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
X <- head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_lt,CHcoHmodel16$results$padj),])
X <- rbind(X,head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_gt,CHcoHmodel16$results$padj),]))

