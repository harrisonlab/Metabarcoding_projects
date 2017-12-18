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

# Convert sequencing run date to a factor
colData$run<-as.factor(colData$run)

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000) & (colData$area=="Rhizosphere"|colData$area=="Stool bed")
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
# plot cummulative reads (will also produce a data table "dtt" in the global environment)
plotCummulativeReads(counts(dds,normalize=T),plot=F)

#### Select filter ####
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

### remove sequencing run bias and plot
pc.res <- resid(aov(mypca$x~colData$run,colData))
df <- t(data.frame(t(pc.res*mypca$percentVar)))

# plot the PCA
ggsave(paste(RHB,"PCA.pdf",sep="_"),plotOrd(df,colData,design="genotype",shape="area",xlabel="PC1",ylabel="PC2",ylims=c(-1.5,0.5)))

#===============================================================================
#       ANOVA
#===============================================================================

# Want to see if genotype can describe more of the variance in the rhizosphere data compared to the bulk soil data.
# This would suggest that the rhizospere are recruting specific OTUs

qfun <- function (mypca,dds,names=c("area","genotype","interaction","residual"),model) {
	X <- t(apply(mypca$x,2,function(x){t(summary(aov(as.formula(paste0("x",model)),colData(dds)))[[1]][2])}))
	colnames(X) <- names
	x<-t(apply(X,1,prop.table))
	perVar <- x * mypca$percentVar
	print(colSums(perVar))
	print(colSums(perVar)/sum(colSums(perVar))*100)
	return(as.data.table(perVar))
}

#pc.res <- resid(aov(mypca$x~run,colData(dds)))

#ss_all_res <- qfun(list(x=pc.res,percentVar=mypca$percentVar),dds)

ss_all <- qfun(mypca,dds,model="~area + genotype + area * genotype",names=c("area","genotype","interaction","residual"))
ss_all_run <- qfun(mypca,dds,model="~run+area + genotype + area * genotype",names=c("run","area","genotype","interaction","residual"))
lapply(seq(1:4),function(x) {X<-summary(aov(mypca$x[,x]~area + genotype + area * genotype,colData(dds)))[[1]][[2]];X/sum(X)*100})
lapply(seq(1:4),function(x) {X<-summary(aov(mypca$x[,x]~run+area + genotype + area * genotype,colData(dds)))[[1]][[2]];X/sum(X)*100})

dds_rhiz  <- dds[,dds$area=="Rhizosphere"]
dds_stool <- dds[,dds$area!="Rhizosphere"]

mypca_rhiz <- des_to_pca(dds_rhiz)
mypca_stool <- des_to_pca(dds_stool)

# sum of squares
ss_pervar_rhiz  <- qfun(mypca_rhiz,dds_rhiz,names=c("run","genotype","residual"),model="~run+genotype")
ss_pervar_rhiz  <- qfun(mypca_rhiz,dds_rhiz,names=c("genotype","residual"),model="~genotype")
ss_pervar_stool <- qfun(mypca_stool,dds_stool,names=c("run","genotype","residual"),model="~run+genotype")
ss_pervar_stool <- qfun(mypca_stool,dds_stool,names=c("genotype","residual"),model="~genotype")

lapply(seq(1:4),function(x) {X<-summary(aov(mypca_rhiz$x[,x]~run+genotype,colData(dds_rhiz)))[[1]][[2]];X/sum(X)*100})
lapply(seq(1:4),function(x) {X<-summary(aov(mypca_rhiz$x[,x]~genotype,colData(dds_rhiz)))[[1]][[2]];X/sum(X)*100})

lapply(seq(1:4),function(x) {X<-summary(aov(mypca_stool$x[,x]~run+genotype,colData(dds_stool)))[[1]][[2]];X/sum(X)*100})
lapply(seq(1:4),function(x) {X<-summary(aov(mypca_stool$x[,x]~genotype,colData(dds_stool)))[[1]][[2]];X/sum(X)*100})

# alternative - combine both P26 sand and clay into single genotype
levels(dds$genotype)[2:3] <- "M26"
mypca <- des_to_pca(dds)
# then re-run this section
