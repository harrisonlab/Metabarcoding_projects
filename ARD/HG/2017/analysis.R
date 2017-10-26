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

# I've rerun the fungal pipeline, which now uses a different naming convention
colnames(ubiom_FUN$countData) <- gsub("([0-9]?[abcdef])([\\._].*)","\\1",colnames(ubiom_FUN$countData)) # \\2 keeps the second match () group

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

#### combine species at 0.95 (default) confidence (if they are species). Works well for Oomycetes and Fungi

# attach OO objects
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
# list of species with more than one associated OTU
combinedTaxa <- combineTaxa("zOO.taxa")
# show the list
combinedTaxa[,1]
# manual filter list to remove none species (e.g. unknown, Pythium aff)
combinedTaxa <- combinedTaxa[c(-3,-9),]
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

# list of species with more than one associated OTU
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("zNEM.taxa")
combinedTaxa <- combinedTaxa[c(-3,-6),]
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
colData <- colData[names(countData),]
# Depending on how I've produced the files...
# row.names(colData) <- colData$name
# colData <- colData[gsub("\\.","-",sub("_.*","",sub("^X","",names(countData)))),]

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000) & (colData$condition!="C")
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
# Correction from aboslute quantification (BAC only at the moment)
sizeFactors(dds)<-geoMeans(dds)* sapply(colData$copies,function(x) x/mean(colData$copies,na.rm=T))


#===============================================================================
#       Filter data 
#============================================================================

### Pythium specific filter to remove OTUs which are unlikely part of the SAR kingdom
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]

### read accumulation filter
# output pdf file
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))					

#### Select filter ####
# Apply seperately for appropriate data set depending on cut-off chosen from graph
myfilter <- dtt$OTU[1:300] #FUN
myfilter <- dtt$OTU[1:40] # OO
myfilter <- dtt$OTU[1:75] # NEM
myfilter <- dtt$OTU[1:4500]  # BAC

# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

### attempt sequencer run information removal (there are better methods available, but this is not so bad)
colData$run <- as.factor(colData$run)					
pc.res <- resid(aov(mypca$x~colData$run,colData))					
					
# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <- t(data.frame(t(pc.res*mypca$percentVar)))

# plot the PCA
pdf(paste(RHB,"VA.pdf",sep="_"))
plotOrd(df,colData,design="condition",xlabel="PC1",ylabel="PC2")
plotOrd(df,colData,shape="condition",design="distance",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

### remove spatial information (this uses the factor "pair" not the numeric "location") and plot
# convert spatial information to a factor
colData$position <- as.factor(colData$meters)	
pc.res <- resid(aov(mypca$x~run+position,colData))
df <- t(data.frame(t(pc.res2*mypca$percentVar)))

ggsave(paste(RHB,"VA_deloc_v2.pdf",sep="_"),plotOrd(df,colData,shape="condition",design="distance",continuous=T,xlabel="PC1",ylabel="PC2"))					

#===============================================================================
#       differential analysis
#===============================================================================
 
# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1
					
# the full model 
full_design <- ~run + position + condition

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
