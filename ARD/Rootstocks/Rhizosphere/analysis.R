#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
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
colnames(ubiom_BAC$countData) <- gsub("\\.","_",sub("_.*","",colnames(ubiom_BAC$countData)))

# Fungi all in one call
ubiom_FUN <- list(
	countData=read.table("FUN.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 
colnames(ubiom_FUN$countData) <- gsub("\\.","_",sub("_.*","",colnames(ubiom_FUN$countData)))

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

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]
# Depending on how I've produced the files...
# row.names(colData) <- colData$name
# colData <- colData[gsub("\\.","-",sub("_.*","",sub("^X","",names(countData)))),]

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000) & (colData$pair!="C")
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - use geoMeans function if
# every gene contains at least one zero, as cannot compute log geometric means
# sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
sizeFactors(dds) <-geoMeans(dds)
# library(edgeR) # I think anyway
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

#===============================================================================
#       Filter data 
#============================================================================

### read accumulation filter
# output pdf file
pdf(paste(RHB,"OTU_counts.pdf",sep="_"))

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
plotCummulativeReads(counts(dds,normalize=T))

# close pdf
dev.off()

#### Select filter ####
# Apply seperately for appropriate data set depending on cut-off chosen from graph
myfilter <- dtt$OTU[1:500] #FUN
myfilter <- dtt$OTU[1:40] # OO

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
colData$location<-as.number(colData$pair)

# plot the PCA
pdf(paste(RHB,"VA.pdf",sep="_"))
plotOrd(df,colData,design="condition",xlabel="PC1",ylabel="PC2")
plotOrd(df,colData,shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

### remove spatial information (this uses the factor "pair" not the numeric "location") and plot
pc.res <- resid(aov(mypca$x~colData$pair,colData))
df <- t(data.frame(t(pc.res*mypca$percentVar)))

pdf(paste(RHB,"VA_deloc.pdf",sep="_"))
plotOrd(df,colData,shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
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

# contrast (not actually necessary in this case as this would be the default result calculated by results(dds)
# contrast <- c("condition","S","H")
res <- results(dds,alpha=alpha,parallel=T)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff_filtered.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# MA plot
pdf(paste(RHB,"ma_plot.pdf",sep="_"))
plot_ma(res.merge)
dev.off()


#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()