#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(BiocParallel)
library(data.table)
library(plyr)
library(tidyverse)
library(metafuncs)
library(vegan)
library(lmPerm)

#environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)
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
#taxData<-phyloTaxaTidy(taxData,0.65)

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


# Fungi all in one call
ubiom_FUN <- list(
  countData=read.table("FUN.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  RHB="FUN"
)


#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("FUN.taxa",column_order=c(1:7,9:15))
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

colData$parent <- as.factor(colData$parent)

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Sample rarefaction plots
#===============================================================================        

library(grid)
library(gridExtra)
library(viridis)
				    
gfunc <- function(countData,coldata,title) {        
  colData <- colData[names(countData),]

  # remove low count and control samples
  myfilter <- colData$Treatment!="Control"

  # apply filter
  colData <- droplevels(colData[myfilter,])
  countData <- countData[,myfilter]

  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))

  # get cummulative sum of each sample
  DT <- cumsum(DT)    

  # log the count values                            
  DT <- log10(DT)

  # relabel columns
  colnames(DT) <- sub("(X)([0-9]+[HS])(.*)","\\2",colnames(DT))

  # set values larger than maximum for each column to NA
  DT <- data.table(apply(DT,2,function(x) {x[(which.max(x)+1):length(x)]<- NA;x}))
  
  # remove rows with all NA
  DT <- DT[rowSums(is.na(DT)) != ncol(DT), ]
  
  # add a count column to the data table
  DT$x <- seq(1,nrow(DT))
                             
  # melt the data table for easy plotting 
  MDT <- melt(DT,id.vars="x")
			      
  # create an empty ggplot object from the data table
  g <- ggplot(data=MDT,aes(x=x,y=value,colour=variable))

  # remove plot background and etc.
  g <- g + theme_classic_thin() %+replace% theme(legend.position="none",axis.title=element_blank())

  # plot cumulative reads
  g <- g + geom_line(size=1.5) + scale_colour_viridis(discrete=T)

  # add axis lables
  g <- g + ggtitle(title)
  #g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")

  # print the plot
  g
}

invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
g1 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Bacteria")
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
g2 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Fungi")
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
g3 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Oomycetes")
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
g4 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Nematodes")

glegend <- get_legend(g)
ggsave("rarefaction_all.pdf",grid.arrange(g1,g2,g3,g4,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=2))                              


#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# Recreate dds object and don't filter for low counts before running Alpha diversity

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values
ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour=NULL,measures=c("Chao1", "Shannon", "Simpson","Observed")))
ggsave(paste(RHB,"Alpha_Chao1.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Chao1"))) # ,limits=c(0,xxx,"Chao1")
ggsave(paste(RHB,"Alpha_Shannon.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Shannon")))
ggsave(paste(RHB,"Alpha_Simpson.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Simpson")))
ggsave(paste(RHB,"Alpha_Observed.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Observed")))

### permutation based anova on diversity index ranks ###
# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",returnData=T)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,colData,by=c("Samples"="Sample_FB"))) # or sample_on

# perform anova for each index
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
sink()




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
plotOrd(d,colData,design="parent",shape="sex",alpha=0.75,cbPalette=T)
plotOrd(d,colData,design="parent",shape="sex",axes=c(2,3),alpha=0.75,cbPalette=T)














#===============================================================================
#       differential analysis
#===============================================================================

# p value for FDR cutoff
alpha <- 0.1

# the model
design <- ~sex*parent

# add design to dds object
design(dds) <- design

# run model
dds <- DESeq(dds,parallel=F)

# difference between sexes
res <- results(dds,alpha=alpha,parallel=F,contrast=c("sex","M","F"))
		   
# difference between 

# results COD vs Healthy
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","COD","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm1 <- res.merge
       
# results AOD vs Healthy
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","AOD","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"AOD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm2 <- res.merge

# results Remission vs Healthy
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","Remission","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"Rem_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm3 <- res.merge

# results COD vs AOD
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","COD","AOD"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_AOD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm4 <- res.merge

# results COD vs Remission
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","COD","Remission"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_REM_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm5 <- res.merge

# results AOD vs Remission
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","AOD","Remission"))
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
sapply(seq(1,4),function(i) write.table(all.tabs[[i]],paste(RHB,names(res.merge)[i],"diffs.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
       
# chestnuts and bigwood
res <- lapply(list_dds[5:6],results,alpha=alpha)
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
       
# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[ res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa")) 

# Merge

rm <- lapply(rm2,function(rm) rm[,c(1,2,3,7),drop=F])      
suffixes <- paste0("_",names(rm))
test <- merge(merge(rm[[1]],rm[[2]],by="OTU",all=T,suffixes=suffixes[1:2]),
	     merge(rm[[3]],rm[[4]],by="OTU",all=T,suffixes=suffixes[3:4]),
	     all=T )
res.merge <- data.table(inner_join(test,data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,"AOD.diff.txt",sep="\t",quote=F,na="") 
