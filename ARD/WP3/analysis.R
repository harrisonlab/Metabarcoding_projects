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

colnames(A1) <- sub("_.*","",sub("-","\\.",colnames(A1)))
colnames(A2) <- sub("_.*","",sub("-","\\.",colnames(A2)))
colnames(A3) <- sub("_.*","",sub("-","\\.",colnames(A3)))
colnames(A4) <- sub("_.*","",sub("-","\\.",colnames(A4)))

temp <- ubiom_BAC$countData
temp$"#OTU ID" <- rownames(temp)
temp <- as.data.frame(as.data.table(rbind.fill(temp,A3))[,lapply(.SD,sum,na.rm=T),by="#OTU ID"])
rownames(temp) <- temp$"#OTU ID"
temp <- temp[,-1]
ubiom_BAC$countData <- temp

temp <- ubiom_FUN$countData
temp$"#OTU ID" <- rownames(temp)
temp <- as.data.frame(as.data.table(rbind.fill(temp,A1))[,lapply(.SD,sum,na.rm=T),by="#OTU ID"])
rownames(temp) <- temp$"#OTU ID"
temp <- temp[,-1]
ubiom_FUN$countData <- temp 

rm(temp)

# ergh some colnames don't match colData - uppercase rep
# easiest just to convert everything to upper or lower
colnames(ubiom_FUN$countData) <- toupper(colnames(ubiom_FUN$countData))
colnames(ubiom_BAC$countData) <- toupper(colnames(ubiom_BAC$countData))
rownames(ubiom_FUN$colData) <- toupper(rownames(ubiom_FUN$colData))
rownames(ubiom_BAC$colData) <- toupper(rownames(ubiom_BAC$colData))

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

ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN)
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC)

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
