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
library(ggplot2)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
library(Biostrings)
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
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 
rownames(ubiom_FUN$colData) <- ubiom_FUN$colData$name

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
	RHB="NEM"
) 

#===============================================================================
#       Combine species 
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species).

# attach OO objects
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
# list of species with more than one associated OTU
combinedTaxa <- combineTaxa("zOO.taxa")
# show the list
combinedTaxa
# manual filter list to remove none species (e.g. unknown, Pythium aff)
combinedTaxa <- combinedTaxa[c(-3,-4,-7,-9,-10,-11),]
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
# all species in combinedTaxa are combinable
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData

# nematodes
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("zNEM.taxa")
combinedTaxa
combinedTaxa <- combinedTaxa[c(-4,-6,-7),]
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_NEM$countData <- countData
ubiom_NEM$taxData <- taxData

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
# for v2 of FUN data
# colnames(countData) <- paste0(sub(".*_","",names(countData)),"D161020") 
colData <- colData[names(countData),]

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000) & (colData$pair!="Control")
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - using geoMeans function (works better with this data set)
max(geoMeans(dds))/min(geoMeans(dds))
max(sizeFactors(estimateSizeFactors(dds)))/min(sizeFactors(estimateSizeFactors(dds)))
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
sizeFactors(dds) <-geoMeans(dds) 
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

#===============================================================================
#       Collapse replicates
#============================================================================
# the data contains three replicates for each sampling point 
# these are likely to be highly correlated and could mess with the differential analysis
dds <- collapseReplicates2(dds,groupby=paste0(dds$condition,dds$pair),simple=F)
colData <- as.data.frame(colData(dds))

#===============================================================================
#       Filter data 
#============================================================================

# Oomcete specific filter to remove OTUs which are unlikely part of the SAR kingdom
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
# apply filter
dds <- dds[myfilter,]
# nematode
myfilter <- row.names(taxData[as.number(taxData$c_conf)>0.9 & as.number(taxData$o_conf)>0.9,])
dds <- dds[rownames(dds)%in%myfilter,]

### read accumulation filter

# plot cummulative reads
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

# Apply seperately for appropriate data set depending on cut-off chosen from graph
myfilter <- dtt$OTU[1:500] #FUN
myfilter <- dtt$OTU[1:80] # OO
myfilter <- dtt$OTU[1:60] # NEM
myfilter <- dtt$OTU[1:5000]  # BAC

# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot 
dds$location<-as.number(dds$pair)

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
plotOrd(df,colData(dds),design="condition",xlabel="PC1",ylabel="PC2")
plotOrd(df,colData(dds),shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

### remove spatial information (this uses the factor "pair" not the numeric "location") and plot
pc.res <- resid(aov(mypca$x~pair,colData(dds)))
df <- t(data.frame(t(pc.res*mypca$percentVar)))

pdf(paste(RHB,"PCA_deloc.pdf",sep="_"))
plotOrd(df,colData(dds),shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================
 
# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

# the full model 
full_design <- ~pair + condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# get results (no contrast needed as default is correct way round)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy table
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# output results
write.table(res.merge, paste(RHB,"diff_filtered.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
writeXStringSet(readDNAStringSet(paste0(RHB,".zotus.fa"))[res.merge[padj<=0.1]$OTU],paste(RHB,".sig.fa"))
