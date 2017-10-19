#===============================================================================
#       Load libraries
#===============================================================================

# library(phyloseq)
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
# library(cooccur)
# library(parallel)

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
colnames(ubiom_FUN$countData) <- gsub("(^.*_)(S[0-9]*)($)","\\2D161020",colnames(ubiom_FUN$countData)) # \\2 keeps the second match () group

# Oomycetes
ubiom_OO <- list(
	countData=read.table("OO.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData2",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zOO.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="OO"
) 

# Nematodes 
ubiom_NEM <- list(
	countData=read.table("NEM.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData2",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("zNEM.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="NEM"
) 

#===============================================================================
#       Combine species 
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species). Works well for Oomycetes and Fungi

# attach OO objects
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
# list of species with more than one associated OTU
combinedTaxa <- combineTaxa("zOO.taxa")
# show the list
combinedTaxa[,1]
# manual filter list to remove none species (e.g. unknown, Pythium aff)
combinedTaxa <- combinedTaxa[c(-3,-9),]
# adjust countData for combined taxa
countData <- combCounts(combinedTaxa,countData)
# adjust taxData for combined taxa
taxData <- combTaxa(combinedTaxa,taxData)
# copy back to object
ubiom_OO$countData <- countData
ubiom_OO$taxData <- taxData

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
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# ensure colData rows and countData columns have the same order
# for v2 of FUN data
# colnames(countData) <- paste0(sub(".*_","",names(countData)),"D161020") 
colData <- colData[names(countData),]

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000) & (colData$pair!="Control")
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# collapse replicates to their mean (collapseReplicates calulates the sum of replicates) - probably best to do this before library size correction
dds$group <- paste(dds$condition,dds$pair,sep="_")
dds <- collapseReplicates2(dds,groupby=dds$group,simple=T)

# calculate size factors - using geoMeans function (works better with this data set)
max(geoMeans(dds))/min(geoMeans(dds))
max(sizeFactors(estimateSizeFactors(dds)))/min(sizeFactors(estimateSizeFactors(dds)))
# sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
sizeFactors(dds) <-geoMeans(dds) 
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

#===============================================================================
#       Filter data 
#============================================================================

# Oomcete specific filter to remove OTUs which are unlikely part of the SAR kingdom
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
# apply filter
dds <- dds[myfilter,]

### read accumulation filter
# output pdf file
pdf(paste(RHB,"OTU_counts.pdf",sep="_"))

# plot cummulative reads
plotCummulativeReads(counts(dds,normalize=T))

# close pdf
dev.off()

#### Select filter ####
# get row sum of normalized counts
df <- as.data.table(rowSums(counts(dds,normalize=T)),keep.rownames=T)
# order decending
df <- df[order(-V2)]

# Apply seperately for appropriate data set depending on cut-off chosen from graph
myfilter <- df$V1[1:700] #FUN
myfilter <- df$V1[1:50] # OO
myfilter <- df$V1[1:80] # NEM
myfilter <- df$V1[1:600]  # BAC

# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot 
dds$location<-as.number(dds$pair)

# plot the PCA
pdf(paste(RHB,"VA.pdf",sep="_"))
plotOrd(df,colData(dds),design="condition",xlabel="PC1",ylabel="PC2")
plotOrd(df,colData(dds),shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

### remove spatial information (this uses the factor "pair" not the numeric "location") and plot
pc.res <- resid(aov(mypca$x~pair,colData(dds)))
df <- t(data.frame(t(pc.res*mypca$percentVar)))

pdf(paste(RHB,"VA_deloc.pdf",sep="_"))
plotOrd(df,colData(dds),shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================
 
# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

# the full model 
full_design <- ~pair + condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# get results (no contrast needed as default is correct way round)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy table
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# output resules
write.table(res.merge, paste(RHB,"diff_filtered.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# MA plot
pdf(paste(RHB,"ma_plot.pdf",sep="_"))
plot_ma(res.merge)
dev.off()

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()

#===============================================================================
#       network analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

cotables <- lapply(myfiltbioms,function(obj) as.data.frame(as.matrix(otu_table(obj))))
cotables_h <- lapply(seq(1,length(myfiltbioms)),function(i)  cotables[[i]][,row.names(sample_data(prune_samples(sample_data(myfiltbioms[[i]])$condition=="Healthy",myfiltbioms[[i]])))]) 
cotables_s <- lapply(seq(1,length(myfiltbioms)),function(i)  cotables[[i]][,row.names(sample_data(prune_samples(sample_data(myfiltbioms[[i]])$condition=="Symptom",myfiltbioms[[i]])))])

cotables_h <- lapply(cotables_h, function(obj) obj[rowSums(obj)>5,colSums(obj)>5])
cotables_s <- lapply(cotables_s, function(obj) obj[rowSums(obj)>5,colSums(obj)>5])

lapply(seq(1,length(cotables_h)), function(i) cotables_h[[i]][cotables_h[[i]]>0] <<- 1)
lapply(seq(1,length(cotables_s)), function(i) cotables_s[[i]][cotables_s[[i]]>0] <<- 1)

CFcoHmodels <- mclapply(cotables_h, function(obj) cooccur2(obj,type = "spp_site",spp_names = T,thresh = T),mc.cores=4)
CFcoSmodels <- mclapply(cotables_s, function(obj) cooccur2(obj,type = "spp_site",spp_names = T,thresh = T),mc.cores=4)

lapply(seq(1,length(CFcoHmodels)), function(i) {
	CFcoHmodels[[i]]$results$padj <<- p.adjust(apply(CFcoHmodels[[i]]$results[,8:9],1, min),"BH")
	CFcoSmodels[[i]]$results$padj <<- p.adjust(apply(CFcoSmodels[[i]]$results[,8:9],1, min),"BH")
})

lapply(CFcoHmodels, function(obj) {nrow(obj$results[obj$results$padj<=0.1,]})

lapply(CFcoHmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_gt<=0.05,]))
lapply(CFcoHmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_lt<=0.05,]))

lapply(CFcoSmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_gt<=0.05,]))
lapply(CFcoSmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_lt<=0.05,]))


X <- rbind.fill(lapply(CFcoHmodels, function(obj) head(obj$results[order(obj$results$padj),],6)))
Y <- rbind.fill(lapply(CFcoSmodels, function(obj) head(obj$results[order(obj$results$padj),],6)))


HcoHmodel16$results$p_lt

write.table(X,"cooc.healthy.txt",sep="\t",quote=F,row.names=F)

head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
X <- head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_lt,CHcoHmodel16$results$padj),])
X <- rbind(X,head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_gt,CHcoHmodel16$results$padj),]))

