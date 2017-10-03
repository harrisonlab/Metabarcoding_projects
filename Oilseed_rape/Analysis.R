#===============================================================================
#       Load libraries 
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
library(Biostrings)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")

#===============================================================================
#       Load data 
#===============================================================================

# load denoised otu count table
countData <- read.table("16S.zotus_table.txt",header=T,sep="\t",row.names=1, comment.char = "")

# load sample metadata
colData <- read.table("colData",header=T,sep="\t",row.names=1)

# load taxonomy data
taxData <- read.table("z16S.taxa",header=F,sep=",",row.names=1)

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
	countData=read.table("ITS.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zITS.taxa",header=T,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# attach objects (FUN or BAC)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# The seqeuncing run contains additional data, so subset rape data
dds <- dds[,dds$type=="rape"]

# remove low count samples
filter <- colSums(counts(dds))>=1000
dds <- dds[,filter]

# calculate size factors - use geoMeans function if
# every gene contains at least one zero, as cannot compute log geometric means
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
#sizeFactors(dds) <-geoMeans(dds)
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
pdf(paste(RHB,"rape.pdf",sep="."))
plotOrd(df,dds@colData,design="condition")
dev.off()

#===============================================================================
#       Differential analysis of rape data
#===============================================================================

# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# drop unused levels from condition
dds$condition <- droplevels(dds$condition)
dds$farm <- droplevels(dds$farm)
dds$field <- droplevels(dds$field)
dds$bean <- droplevels(dds$bean)

# p value for FDR cutoff
alpha <- 0.1

# add model to the DES object
design(dds) <- ~condition

# calculate fit
dds <- DESeq(dds,parallel=T)

# main effect
contrast=c("condition","low_till", "cultivated" )

# calculate results
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast,cooksCutoff=F)

# merge results with taxonomy table
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# show those with BH p <= 0.05
res.merge[padj<=0.05,]

write.table(res.merge,paste(RHB,"rape_main_effect.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

