#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(BiocParallel)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(Biostrings)
library(vegan)
library(lmPerm)
library(phyloseq)
library(ape)

register(MulticoreParam(12))
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)
#assignInNamespace("plot_ordination",value=plot_ordination,ns="phyloseq")

#===============================================================================
#       Load data
#===============================================================================

# load otu count table
countData <- read.table("BAC.otus_table.txt",header=T,sep="\t",row.names=1, comment.char = "")

# load sample metadata
colData <- read.table("colData",header=T,sep="\t",row.names=1)

# load taxonomy data
taxData <- read.table("BAC.taxa",header=F,sep=",",row.names=1)

# reorder columns
taxData<-taxData[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]

# add best "rank" at 0.65 confidence and tidy-up the table
taxData<-phyloTaxaTidy(taxData,0.65)

# get unifrac dist
phylipData <- fread.phylip("BAC.phy")

njtree <- nj(as.dist(phylipData))

# save data into a list
ubiom_BAC <- list(
  countData=countData,
  colData=colData,
  taxData=taxData,
  phylipData=phylipData,
  njtree=njtree,
  RHB="BAC"
)

# Fungi all in one call
ubiom_FUN <- list(
  countData=read.table("FUN.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  phylipData=fread.phylip("FUN.phy"),
  RHB="FUN"
)
ubiom_FUN$njtree <- nj(as.dist(ubiom_FUN$phylipData))

# Oomycetes
ubiom_OO <- list(
  countData=read.table("OO.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData2",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("OO.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  phylipData=fread.phylip("OO.phy"),
  RHB="OO"
)
ubiom_OO$njtree <- nj(as.dist(ubiom_OO$phylipData))
rownames(ubiom_OO$colData) <- paste0("X",gsub("_","\\.",ubiom_OO$colData$name),"_",sub("D.*","",rownames(ubiom_OO$colData)))
rownames(ubiom_OO$colData) <- sub("XG","G",rownames(ubiom_OO$colData))

# Nematodes
ubiom_NEM <- list(
  countData=read.table("NEM.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData2",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("NEM.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  phylipData=fread.phylip("NEM.phy"),
  RHB="NEM"
)
ubiom_NEM$njtree <- nj(as.dist(ubiom_NEM$phylipData))
rownames(ubiom_NEM$colData) <- paste0("X",gsub("_","\\.",ubiom_NEM$colData$name),"_",sub("D.*","",rownames(ubiom_NEM$colData)))
rownames(ubiom_NEM$colData) <- sub("XG","G",rownames(ubiom_NEM$colData))

#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("FUN.taxa")
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData

# Nematodes

# oomycetes
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
#rownames(colData) <- sub("^XG","G",rownames(colData))
colData <- colData[names(countData),]

# remove low count and control samples 
filter <- (colSums(countData)>=1000) & colData$condition!="C"

# remove pair of any sample with a low count
exclude<-which(!filter)
filter <- filter&sapply(colData$pair,function(x) length(which(x==colData$pair[-exclude]))>1)

# apply filter
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - use geoMeans function if every gene contains at least one zero, as cannot compute log geometric means (or calcNormFactors)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)
# library(edgeR) 
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))
