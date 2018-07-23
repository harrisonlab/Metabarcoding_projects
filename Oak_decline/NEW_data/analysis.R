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
#phylipData <- fread.phylip("BAC.phy")

#njtree <- nj(as.dist(phylipData))

# save data into a list
ubiom_BAC <- list(
  countData=countData,
  colData=colData,
  taxData=taxData,
 # phylipData=phylipData,
  #njtree=njtree,
  RHB="BAC"
)
rownames(ubiom_BAC$colData) <- sub(".*_L1_","",rownames(ubiom_BAC$colData))
names(ubiom_BAC$countData) <- sub(".*_L1_","",names(ubiom_BAC$countData))

# Fungi all in one call
ubiom_FUN <- list(
  countData=read.table("FUN.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
#  phylipData=fread.phylip("FUN.phy"),
  RHB="FUN"
)
rownames(ubiom_FUN$colData) <- sub(".*_L1_","",rownames(ubiom_FUN$colData))
names(ubiom_FUN$countData) <- sub(".*_L1_","",names(ubiom_FUN$countData))

#ubiom_FUN$njtree <- nj(as.dist(ubiom_FUN$phylipData))

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
#      ****FUNGI/BACTERIA****
#===============================================================================

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
# Bacteria
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
# ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================

# p value for FDR cutoff
alpha <- 0.1

# the model
design <- ~condition

# split dds object into per wood
# first get rid of bigwood as it only has 2 samples - and remove it from the levels of site
list_dds <-list(Attingham   = dds[,dds$site=="Attingham"],
		            Gt_Monk     = dds[,dds$site=="Gt_Monk"],
		            Langdale    = dds[,dds$site=="Langdale"],
		            Winding     = dds[,dds$site=="Winding"],
		            Chestnuts   = dds[,dds$site=="Chestnuts"],
		            Chestnuts   = dds[,dds$site=="Bigwood"])		
# Filter low count OTUs
list_dds <- lapply(list_dds,function(dds)  dds[rowSums(counts(dds, normalize=T))>0,])

# add full model to dds object
list_dds <- lapply(list_dds,function(dds) {
	design(dds) <- design;
	colData(dds) <- droplevels(colData(dds))
	dds}
)

# calculate fit
list_dds <- lapply(list_dds,DESeq,parallel=T)

# results COD vs AOD
res <- lapply(list_dds[1:4],results,alpha=alpha,parallel=T,contrast=c("condition","COD","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm1 <- res.merge
       
# results COD vs AOD
res <- lapply(list_dds[1:4],results,alpha=alpha,parallel=T,contrast=c("condition","AOD","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"AOD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm2 <- res.merge

# results COD vs AOD
res <- lapply(list_dds[1:4],results,alpha=alpha,parallel=T,contrast=c("condition","Remission","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"Rem_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm3 <- res.merge

# results COD vs AOD
res <- lapply(list_dds[1:4],results,alpha=alpha,parallel=T,contrast=c("condition","COD","AOD"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_AOD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm4 <- res.merge

# results COD vs Remission
res <- lapply(list_dds[1:4],results,alpha=alpha,parallel=T,contrast=c("condition","COD","Remission"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_REM_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm5 <- res.merge

# results AOD vs Remission
res <- lapply(list_dds[1:4],results,alpha=alpha,parallel=T,contrast=c("condition","AOD","Remission"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"AOD_REM_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm6 <- res.merge
       
all.tabs <- lapply(seq(1,4),function(i) data.table(list(
	rm1[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_COD = log2FoldChange)  %>% rename(padj_COD = padj),
	rm2[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_AOD = log2FoldChange)  %>% rename(padj_AOD = padj),
	rm3[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_REM = log2FoldChange)  %>% rename(padj_REM = padj),
	rm4[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_CA  = log2FoldChange)  %>% rename(padj_CA  = padj),
	rm5[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_CR  = log2FoldChange)  %>% rename(padj_CR  = padj),
	rm6[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_AR  = log2FoldChange)  %>% rename(padj_AR  = padj)
        ) %>% Reduce(function(dt1,dt2) full_join(dt1,dt2), .))
)	
sapply(seq(1,4),function(x) write.table(all.tabs[[x]],paste(RHB,names(res.merge)[x],"diffs.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
       
# chestnuts and bigwood
res <- lapply(list_dds[5:6],results,alpha=alpha)
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
       
# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[ res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa")) 
