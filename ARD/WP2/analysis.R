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


ubiom_BAC <- loadData("BAC.otus_table.txt","colData","BAC.taxa","BAC.phy",RHB="BAC")
ubiom_FUN <- loadData("FUN.otus_table.txt","colData","FUN.taxa","FUN.phy",RHB="FUN")
ubiom_OO <- loadData("OO.otus_table.txt","colData","OO.taxa","OO.phy",RHB="OO")
ubiom_NEM <- loadData("NEM.otus_table.txt","colData","NEM.taxa","NEM.phy",RHB="NEM")


ubiom_BAC$njtree <- nj(as.dist(ubiom_BAC$phylipData))
ubiom_FUN$njtree <- nj(as.dist(ubiom_FUN$phylipData))
ubiom_OO$njtree <- nj(as.dist(ubiom_OO$phylipData))
ubiom_NEM$njtree <- nj(as.dist(ubiom_NEM$phylipData))


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
