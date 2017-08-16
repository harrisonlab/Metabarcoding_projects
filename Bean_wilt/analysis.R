#library(phyloseq)
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(dplyr)
library(plyr)
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

# save data into a list, then ubiom_16S$countData to access countData and etc.
ubiom_16S <- list(
	countData=countData,
	colData=colData[colData,],
	taxData=taxData
)

# remove 18S level from colData
ubiom_16S$colData$Loci <- droplevels(ubiom_16S$colData$Loci)

# or all in one
ubiom_ITS <- list(
	countData=read.table("FUN.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=colData[colData,],
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=T,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65)
) 

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# attach either ITS or 16S objects
invisible(mapply(assign, names(ubiom_ITS), ubiom_ITS, MoreArgs=list(envir = globalenv())))


# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - using geoMeans function due to DES error
# every gene contains at least one zero, cannot compute log geometric means
#sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
sizeFactors(dds) <-geoMeans(dds)

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decompossion on DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
pdf("all.its.pdf")
plotOrd(df,dds@colData,design="condition",shape="type")
dev.off()

# bean sample X45_S68 is a clear outlier

#===============================================================================
#       Differential analysis of bean data
#===============================================================================

# subset bean data
dds <- dds[,dds$type=="bean"]

# remove ITS outlier as identified by PCA (no 16S outliers)
dds <- dds[,colnames(dds)!="X45_S68"]


# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# dds<-dds[rowSums(counts(dds,normalize=T))>0,]


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

# contrasts to test
# main effect
contrast=c("condition","SICK", "HEALTHY" )
contrast=c("bean","FRENCH","RUNNER")

# interaction term (both should give roughly the same results - with FC in opposite direction)
#contrast=list( "beanFRENCH.conditionHEALTHY","beanFrench.conditionHEALTHY") # effect of condition on french beans 
contrast=list( "beanRUNNER.conditionHEALTHY","beanRUNNER.conditionSICK") # effect of condition on runner beans

# calculate results
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)

# merge results with taxonomy table
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# show those with BH p <= 0.05
res.merge[padj<=0.05,]

write.table(res.merge,"main_effect.txt",quote=F,sep="\t",na="",row.names=F)
write.table(res.merge,"interaction.txt",quote=F,sep="\t",na="")
write.table(res.merge,"bean_effect.txt",quote=F,sep="\t",na="")

#===============================================================================
#       Differential OTUs
#===============================================================================

# read OTU fasta
myOTUs <- readDNAStringSet("ITS.otus.fa")

writeXStringSet(myOTUs[as.matrix(res.merge[padj<=0.1,1]),"main_big.fa")

#===============================================================================
#       PCA plots of beans only
#===============================================================================

mypca <- des_to_pca(dds)
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

plotOrd(df,dds@colData,design="condition",title="Condition")
plotOrd(df,dds@colData,design="farm",title="Farm")
plotOrd(df,dds@colData,design="variety", shape="bean",title="Bean Type")
plotOrd(df,dds@colData,design="field", shape="variety",title="Field")



