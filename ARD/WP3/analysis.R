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

#===============================================================================
#       Load data
#===============================================================================

ubiom_BAC <- loadData("BAC.otu_table.txt","colData","BAC.taxa",RHB="BAC")
ubiom_FUN <- loadData("FUN.otu_table.txt","colData","FUN.taxa",RHB="FUN")

A1 <- fread("ambiguous1.otu_table.txt") # fungi r1
A2 <- fread("ambiguous2.otu_table.txt") # bacteria r1 (not used at moment)
A3 <- fread("ambiguous3.otu_table.txt") # bacteria merged
A4 <- fread("ambiguous4.otu_table.txt") # fungi merged (not used)

ubiom_BAC$countData <- as.data.table(rbind.fill(ubiom_BAC$countData,A3))[,lapply(.SD,sum,na.rm=T),by="#OTU ID"] 
ubiom_FUN$countData <- as.data.table(rbind.fill(ubiom_BAC$countData,A1))[,lapply(.SD,sum,na.rm=T),by="#OTU ID"] 

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

#===============================================================================
#       Create DEseq objects
#===============================================================================

ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN,filter=expression(colSums(countData)>=1000&colData$Block!="R"))
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC,filter=expression(colSums(countData)>=1000&colData$Block!="R"))

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
