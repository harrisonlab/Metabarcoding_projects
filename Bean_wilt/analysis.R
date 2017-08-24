#===============================================================================
#       Load libraries 
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(dplyr)
library(Biostrings)
library(devtools)
load_all("../metabarcoding_pipeline/scripts/myfunctions")

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
	colData=colData[colData,],
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

# Oomycetes all in one call
ubiom_OO <- list(
	countData=read.table("OO.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colDataOO",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zOO.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="OO"
) 

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# attach objects (FUN, BAC or OO)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
#invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
#invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# The seqeuncing run contains additional data, so subset bean data
dds <- dds[,dds$type=="bean"]

# calculate size factors - use geoMeans function if
# every gene contains at least one zero, as cannot compute log geometric means
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
#sizeFactors(dds) <-geoMeans(dds)

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
pdf(paste(RHb,"pdf",sep="."))
plotOrd(df,dds@colData,design="condition",)
dev.off()

# Fungi sample X45_S68 is a clear outlier 

#===============================================================================
#       Differential analysis of bean data
#===============================================================================

# remove FUN outlier as identified by PCA (no BAC or OO outliers)
dds <- dds[,colnames(dds)!="X45_S68"]

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
design(dds) <- ~farm + bean + condition + condition:bean

# calculate fit
dds <- DESeq(dds,parallel=T)

# main effect
contrast=c("condition","SICK", "HEALTHY" )

# bean effect
contrast=c("bean","FRENCH","RUNNER")

# interaction term (both should give roughly the same results - with FC in opposite direction)
# contrast=list( "beanFRENCH.conditionHEALTHY","beanFrench.conditionHEALTHY") # effect of condition on french beans 
contrast=list( "beanRUNNER.conditionHEALTHY","beanRUNNER.conditionSICK") # effect of condition on runner beans

# calculate results
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)

# merge results with taxonomy table
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# show those with BH p <= 0.05
res.merge[padj<=0.05,]

write.table(res.merge,paste(RHB,"main_effect.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
write.table(res.merge,paste(RHB,"interaction.txt",sep="_"),quote=F,sep="\t",na="")
write.table(res.merge,paste(RHB,"bean_effect.txt",sep="_"),quote=F,sep="\t",na="")