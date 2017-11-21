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

#### combine species at 0.95 (default) confidence (if they are species). Works well for Oomycetes and Fungi

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

# attach objects (FUN, BAC)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# remove low count samples and control samples (not needed here)
filter <- colSums(countData)>=1000
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - use geoMeans function if
# every gene contains at least one zero, as cannot compute log geometric means
 sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)
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
myfilter <- dtt$OTU[1:4500]  # BAC

myfilter <- dtt$OTU[dtt$CD>5]
# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
plotOrd(df,colData,design="condition",shape="block",xlabel="PC1",ylabel="PC2")

### remove/minimise block effect
pc.res <- resid(aov(mypca$x~colData$block,colData))
df <- t(data.frame(t(pc.res*mypca$percentVar)))
plotOrd(df,colData,design="condition",shape="block",xlabel="PC1",ylabel="PC2")
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
full_design <- ~block + condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

contrast <- c("condition","Treated","Untreated") # the default  alculates the contrast the othe rway
res <- results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff_filtered.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

#################################################################################
#################################################################################

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()
