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

register(MulticoreParam(12))
load_all("~/pipelines/metabarcoding/scripts/myfunctions")

#===============================================================================
#       Load data and tidy up names
#===============================================================================

ubiom_FUN <- loadData("FUN.otus_table.txt","colData","FUN.taxa",RHB="FUN")
colnames(ubiom_FUN$countData) <- gsub("\\.","_",sub("_.*","",colnames(ubiom_FUN$countData)))
ubiom_OO <- loadData("OO.otus_table.txt","colData","OO.taxa",RHB="OO")
colnames(ubiom_OO$countData) <- gsub("\\.","_",sub("_.*","",colnames(ubiom_OO$countData)))

#===============================================================================
#       Combine species
#===============================================================================
# none to combine

#===============================================================================
#       Create DEseq objects
#===============================================================================

ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN,filter=expression(colSums(countData)>=1000),calcFactors=geoMeans)
ubiom_OO$dds <- ubiom_to_des(ubiom_OO,filter=expression(colSums(countData)>=1000&colData$FARM!="DDT TRIAL SITE"),calcFactors=geoMeans)

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, or OO)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Filter data 
#===============================================================================

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

### filter to remove OTUs which are unlikely part of the correct kingdom (best to do this before Alpha diversity analysis)
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]

#===============================================================================
#       Beta diversity
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
ggsave(paste(RHB,"PCA.pdf",sep="_"),plotOrd(d,colData(dds),shape=c("VARIETY"),design="STATUS",continuous=F,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))

#===============================================================================
#       differential analysis
#===============================================================================

# p value for FDR cutoff
alpha <- 0.1

# the full model
full_design <- ~VARIETY + STATUS

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds  <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T,contrast=c("STATUS","SICK","HEALTHY"))

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"VARIETY_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.1]$OTU],paste0(RHB,".VARIETY.sig.fa"))

## Within bean
dds2 <- dds[,dds$BEAN=="RUNNER"]
colData(dds2) <- droplevels(colData(dds2))
design(dds2) <- ~ FARM + STATUS
dds2 <- DESeq(dds2)

res <- results(dds,alpha=alpha,parallel=T)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"RUNNER_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.1]$OTU],paste0(RHB,".RUNNER.sig.fa"))
