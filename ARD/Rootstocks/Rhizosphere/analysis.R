#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")

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
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 

#===============================================================================
#       Combine species 
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# list of species with more than one associated OTU
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("zFUN.taxa")
# all species in combinedTaxa are combinable
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

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# ensure colData rows and countData columns have the same order, then rename them (removing sample IDs)
colData <- colData[colnames(countData),]
colnames(countData) <- rownames(colData) <- colData$name 

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000)
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - further option methods given
 sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)
# library(edgeR) # I think anyway
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

#===============================================================================
#       Filter data 
#===============================================================================

### read accumulation filter
# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

#### Select filter ####
myfilter <- dtt$OTU[dtt$CD>5]
# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA plot all data
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
ggsave(paste(RHB,"PCA.pdf",sep="_"),plotOrd(df,colData,design="genotype",shape="run",xlabel="PC1",ylabel="PC2"))

### remove run information (can't distinguish between run and site) and plot
pc.res <- resid(aov(mypca$x~colData$run,colData))
df <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_deloc.pdf",sep="_"),plotOrd(df,colData,shape="run",design="genotype",xlabel="PC1",ylabel="PC2"))

#===============================================================================
#       Matched genotypes
#===============================================================================

filter <- colData$genotype=="M9"|like(colData$genotype,"M26")
colData <- droplevels(colData[filter,])
countData <- countData[,filter]
dds <- dds[,filter]
dds$genotype <- droplevels(dds$genotype)
dds$run <- droplevels(dds$run)

# pca plot
mypca <- des_to_pca(dds)
df <-t(data.frame(t(mypca$x)*mypca$percentVar))
pc.res <- resid(aov(mypca$x~colData$run,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
pdf(paste(RHB,"matched.pdf",sep="_"))
 plotOrd(df,colData,design="genotype",shape="run",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,design="genotype",shape="run",xlabel="PC1",ylabel="PC2")
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================
 
# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

# add a condition column (this is actually the genotype and genotype is conditon, but it doesn't matter so much)
dds$condition <- as.factor(sub(" .*","",dds$genotype))

dds$gs <- as.factor(paste(dds$genotype,dds$run,sep="_"))

design <- ~gs

# add full model to dds object
design(dds) <- design

# calculate fit
dds <- DESeq(dds,parallel=T)
# contrast (not actually necessary in this case as this would be the default result calculated by results(dds)
# contrast <- c("condition","S","H")
contrast_list <- list(M9=c("gs","M9_a","M9_b"),M26_a=c("gs","M26_a","M26.Sand_b"),M26_b=c("gs","M26_a","M26.Clay_b") ,M26_a=c("gs","M26.Sand_b","M26.Clay_b"))
res_list <- lapply(contrast_list,function(l) {results(dds,alpha=alpha,parallel=T,contrast=l)})
sapply(res_list,summary)
counter<-1
res.merge_list <- lapply(res_list,function(l) {
	X<-data.table(inner_join(data.table(OTU=rownames(l),as.data.frame(l)),data.table(OTU=rownames(taxData),taxData)))
	write.table(X, paste(RHB,names(res_list)[counter],"diff_filtered.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
	counter<<-counter+1
})
dds2 <-dds
design <- ~run+condition
design(dds2) <- design
dds2 <- DESeq(dds2,parallel=T)
res <-results(dds2,alpha=alpha,parallel=T)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"m9-m27_diff_filtered.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()
