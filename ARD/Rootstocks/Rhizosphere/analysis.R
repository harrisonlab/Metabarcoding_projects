#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
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

# ensure colData rows and countData columns have the same order, then rename them (removing sample IDs)
colData <- colData[colnames(countData),]
colnames(countData) <- rownames(colData) <- colData$name 

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000)
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - further option methods given
 sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)
# library(edgeR) # I think anyway
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

# add a condition (actual genotype) column 
dds$condition <- as.factor(sub(" .*","",dds$genotype))
dds$site <- as.factor(sub("_.*","",dds$name))
#levels(dds$site) <- c("Frank Matthews","Shrama")

#===============================================================================
#       Filter data 
#===============================================================================

### read accumulation filter
# plot cummulative reads (will also produce a data table "dtt" in the global environment)
plotCummulativeReads(counts(dds,normalize=T),plot=F)

#### Select filter ####
myfilter <- dtt$OTU[dtt$CD>5]
# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA plot all data
#===============================================================================

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
ggsave(paste(RHB,"PCA.pdf",sep="_"),plotOrd(df,colData(dds),design="genotype",shape="site",xlabel="PC1",ylabel="PC2"))

### remove run information (can't distinguish between run and site) and plot
pc.res <- resid(aov(mypca$x~colData$run,colData))
df <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_deloc.pdf",sep="_"),plotOrd(df,colData(dds),shape="site",design="genotype",xlabel="PC1",ylabel="PC2",xlim=c(-1,1)))

#===============================================================================
#       Matched genotypes across sites
#===============================================================================

dds2 <- dds[,dds$genotype=="M9"|like(dds$genotype,"M26")]
colData(dds2) <- droplevels(colData(dds2))
levels(dds2$genotype)[1:3] <- "M26"

#### pca plot ####
mypca <- des_to_pca(dds2)
df <-t(data.frame(t(mypca$x)*mypca$percentVar))
pc.res <- resid(aov(mypca$x~run,colData(dds2)))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_matched.pdf",sep="_"),plotOrd(df,colData(dds2),design="genotype",shape="site",xlabel="PC1",ylabel="PC2"))
ggsave(paste(RHB,"PCA_matched_no_site.pdf",sep="_"),plotOrd(d,colData(dds2),design="genotype",shape="site",xlabel="PC1",ylabel="PC2"))
#####

#### differential analysis ####
 
# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
dds2<-dds2[rowSums(counts(dds2,normalize=T))>5,]

# p value for FDR cutoff
alpha <- 0.1

# the model (condition=genotype)
design <- ~site+genotype

# add model to dds object
design(dds2) <- design

# calculate
dds2 <- DESeq(dds2,parallel=T)

# extract results
res <- results(dds2,parallel=T,alpha=alpha) # difference is in relation to M9

# combine with taxonomy
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# write to file
write.table(res.merge, paste(RHB,"m9_m26_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
#####

#===============================================================================
#       Within sites (may be better to do this without pooling data)
#===============================================================================
 
# frank mathews
dds_fm <- dds[,dds$site=="FM"]
colData(dds_fm) <- droplevels(colData(dds_fm))
dds_fm <- dds_fm[rowSums(counts(dds_fm,normalize=T))>5,]
# pca plot
mypca <- des_to_pca(dds_fm)
ggsave(paste(RHB,"PCA_FM.pdf",sep="_"),plotOrd(t(data.frame(t(mypca$x)*mypca$percentVar)),colData(dds_fm),design="genotype",xlabel="PC1",ylabel="PC2",ylims=c(5,-5)))
# differential analysis
design(dds_fm) <- ~genotype
dds_fm <- DESeq(dds_fm,reduced=~1,test="LRT",parallel=T)
res_fm_lrt <- results(dds_fm,parallel=T)
res_fm_m9_vs_m25 <- results(dds_fm,parallel=T,names="genotypeM9_vs_genotypeM25")
res_fm_m26_vs_m25 <- results(dds_fm,parallel=T,names="genotypeM26_vs_genotypeM25")
res_fm_m9_vs_m26 <- results(dds_fm,parallel=T,contrast=c("genotype","M9","M26"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res_fm_lrt),as.data.frame(res_fm_lrt)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"FPM_likelihood_ratio_test.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

#===============================================================================
#       differential analysis
#===============================================================================

#===============================================================================
#       Beta diversity analysis
#===============================================================================


#phylip_data<-fread(paste0(RHB,".z.phy"))
#colnames(phylip_data) <- c("OTU",t(phylip_data[,1]))

phylip_data <- fread.phylip(paste0(RHB,".z.phy"))

nj.tree <- nj(as.dist(phylip_data))
write.tree(nj.tree,paste0(RHB,".tree"))

nj.tree <- root(nj.tree,outgroup=200,resolve.root=T)

# the OTU names in the dist object and count data must match - they won't for fungi due to combining similar OTUs
rownames(dds) <- sub("_.*","",rownames(dds))

unifracs <- GUniFrac(t(counts(dds,normalize=T)),nj.tree,alpha=c(0, 0.5, 1))$unifracs

dw <- unifracs[, , "d_1"] # Weighted UniFrac
du <- unifracs[, , "d_UW"] # Unweighted UniFrac

adonis(as.dist(du)~location+Orchard*Sample,colData,parallel=12,permutations=9999))
adonis(as.dist(dw)~run+genotype,colData,parallel=12,permutations=9999)

g <- plotHeatmap(m,textSize=16)
ggsave("test.pdf",g)

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()
