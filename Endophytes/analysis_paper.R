#===============================================================================
#       Load libraries 
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(plyr)
library(dplyr)
library(Biostrings)
library(locfit)
library(ggplot2)
library(viridis)
library(devtools)
library(ape)
library(vegan)
library(phyloseq)

register(MulticoreParam(12))
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)
#
#===============================================================================
#       Load data 
#===============================================================================

# load denoised otu count table
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
# create tree - set.seed to something
sum(utf8ToInt("Greg Deakin"))
njtree <- nj(as.dist(phylipData))
# save data into a list, then ubiom_16S$countData to access countData and etc.
ubiom_BAC <- list(
  countData=countData,
  colData=colData[colData$Loci!="18S",],
  taxData=taxData,
  phylipData=phylipData,
  njtree=njtree,
  RHB="BAC"
)
# remove 18S level from colData
ubiom_BAC$colData$Loci <- droplevels(ubiom_BAC$colData$Loci)

# or all in one
ubiom_FUN <- list(
	countData=read.table("FUN.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=colData[colData$Loci!="16S",],
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	phylipData=fread.phylip("FUN.phy"),
	RHB="FUN"
) 
ubiom_FUN$colData$Loci <- droplevels(ubiom_FUN$colData$Loci)
ubiom_FUN$njtree <- nj(as.dist(ubiom_FUN$phylipData))

#===============================================================================
#       Create DESeq objects 
#===============================================================================

# "attach" the required verion of countData, colData and TaxData
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# set the year column to a factor or deseq won't give expected (correct) results
colData$Year <- as.factor(colData$Year)

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# filter data
dds <- dds[,(colSums(counts(dds))>=1000)&(colData(dds)$Treatment!="Yeast X 2")]

dds$time <- as.integer(sub(" week","",dds$Time.point))
# calculate size factors - use geoMeans function if every gene contains at least one zero (check for size factor range as well)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
#sizeFactors(dds) <-geoMeans(dds) # vst only

# Test size factors

X1 <- colSums(counts(dds))
X2 <- sizeFactors(estimateSizeFactors(dds))
X3 <- geoMeans(dds)
max(X1)/min(X1)
max(X2)/min(X2)
max(X3)/min(X3)

# remove batch effect with limma removeBatchEffect method (maybe useful for plotting 2016/2017 data on same graph)
#library(limma)
#debatched <- removeBatchEffect(counts(dds,normalize=T),colData$Year)

# remove batch effect vst/pca/resid(aov) 
#debatched2 <- batchEffectDes(dds,"Year")


#===============================================================================
#       PCA plot
#===============================================================================

# perform PC decompossion on DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df  <- t(data.frame(t(mypca$x)*mypca$percentVar))

# exclude variability "explained" by Year and Country
pc.res <- resid(aov(mypca$x~Year+Country+Year*Country,dds@colData))

# as above for residual values
d <- t(data.frame(t(pc.res)*mypca$percentVar))

# plot the PCA
pdf(paste(RHB,"saprophyte_v2.pdf",sep="_"))
plotOrd(df,dds@colData,design="Year",shape="Country")
plotOrd(d,dds@colData,design="Year",shape="Country")
plotOrd(df,dds@colData,shape="Treatment",design="Time.point")
plotOrd(d,dds@colData,shape="Treatment",design="Time.point")
dev.off()

# Fungal saprophyte sample X205_S37 doesn't cluster as well as the other samples

# NMDS

# set seed
sum(utf8ToInt("Greg Deakin"))

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object (uses random)
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

#plot 
ggsave(paste(RHB,"Unifrac_NMDS_county_year.pdf",sep="_"),plotOrd(ordu$points,colData(dds),shape="Country",design="Year",xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))
ggsave(paste(RHB,"Unifrac_NMDS_Treatment_time.pdf",sep="_"),plotOrd(ordu$points,colData(dds),shape="Treatment",design="Time.point",xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

#===============================================================================
#       Differential analysis of saprophyte data
#===============================================================================

# drop the Yeast2 level from the Treatment factor (or all unused levels)
colData(dds) <- droplevels(colData(dds))

# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

### Overall effects ####

# the full model 
full_design <- ~Year + Country  + Treatment #+ Time.point + Treatment:Time.point

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# urea
res_urea <- results(dds,contrast=c("Treatment","Urea","Control"),alpha=alpha,parallel=T)
fwrite(data.table(inner_join(data.table(OTU=rownames(res_urea),as.data.frame(res_urea)),data.table(OTU=rownames(taxData),taxData))),paste(RHB,"Urea_overall_effect.txt",sep="_"),sep="\t")

# yeast
res_yeast <- results(dds,contrast=c("Treatment","Yeast","Control"),alpha=alpha,parallel=T)
fwrite(data.table(inner_join(data.table(OTU=rownames(res_yeast),as.data.frame(res_yeast)),data.table(OTU=rownames(taxData),taxData))),paste(RHB,"Yeast_overall_effect.txt",sep="_"),sep="\t")

# urea_yeast
res_urea_yest <- results(dds,contrast=c("Treatment","Urea","Yeast"),alpha=alpha,parallel=T)
fwrite(data.table(inner_join(data.table(OTU=rownames(res_urea_yest),as.data.frame(res_urea_yest)),data.table(OTU=rownames(taxData),taxData))),paste(RHB,"Urea_vs_Yeast_effect.txt",sep="_"),sep="\t")


## difference over time ##

dds_urea       <- dds[,colData(dds)$Treatment!="Yeast"]
dds_yeast      <- dds[,(colData(dds)$Treatment!="Urea")&(colData(dds)$Year!="2016")]
dds_urea_yeast <- dds[,colData(dds)$Treatment!="Control"&(colData(dds)$Year!="2016")]

colData(dds_urea)       <- droplevels(colData(dds_urea))
colData(dds_yeast)      <- droplevels(colData(dds_yeast))
colData(dds_urea_yeast) <- droplevels(colData(dds_urea_yeast))

# add full design to model
design(dds_urea) <- full_design_1 <- ~Year + Country  + Treatment + Time.point + Treatment:Time.point
design(dds_urea_yeast) <- design(dds_yeast) <-  full_design_2 <- ~Country  + Treatment + Time.point + Treatment:Time.point

# the reduced model (for calculating response to time)
reduced_design_1 <- ~Year + Country + Treatment + Time.point
reduced_design_2 <- ~Country + Treatment + Time.point

# calculate model, including both full and reduced designs
dds_urea       <-DESeq(dds_urea, betaPrior=FALSE, test="LRT",full=full_design_1,reduced=reduced_design_1,parallel=T)
dds_yeast      <-DESeq(dds_yeast, betaPrior=FALSE, test="LRT",full=full_design_2,reduced=reduced_design_2,parallel=T)
dds_urea_yeast <-DESeq(dds_urea_yeast, betaPrior=FALSE, test="LRT",full=full_design_2,reduced=reduced_design_2,parallel=T)

# calculate OTUs which respond differently over time (can ignore LFC in output as it's meaningless)	
res_urea_time       <- results(dds_urea,alpha=alpha,parallel=T)	
res_yeast_time      <- results(dds_yeast,alpha=alpha,parallel=T)
res_urea_yeast_time <- results(dds_urea_yeast,alpha=alpha,parallel=T)

fwrite(data.table(inner_join(data.table(OTU=rownames(res_urea_time),as.data.frame(res_urea_time)),data.table(OTU=rownames(taxData),taxData))),paste(RHB,"Urea_time_interaction.txt",sep="_"),sep="\t")
fwrite(data.table(inner_join(data.table(OTU=rownames(res_yeast_time),as.data.frame(res_yeast_time)),data.table(OTU=rownames(taxData),taxData))),paste(RHB,"Yeast_time_interaction.txt",sep="_"),sep="\t")
fwrite(data.table(inner_join(data.table(OTU=rownames(res_urea_yeast_time),as.data.frame(res_urea_yeast_time)),data.table(OTU=rownames(taxData),taxData))),paste(RHB,"Urea_vs_Yeast_time_interaction.txt",sep="_"),sep="\t")


#===============================================================================
#       Graph analysis
#===============================================================================

# filter out rows with mean less than 1 (could probably go higher)
dds2 <-dds[rowMeans(counts(dds,normalize=T))>=1,]
		    
#design(dds2) <- ~farm + date + condition

# calculate rld (log value)
rld <- rlog(dds2,blind=F)

# order results by largest row sum first
rld <- rld[order(rowSums(assay(rld)),decreasing=T),]

# or use vst - not so good if size factors differ markedly	    
vst <- varianceStabilizingTransformation(dds2)		    

X<-unique(c(AS,US,YS))
vst <- vst[X[X%in%row.names(dds2)],]		    
vst <- vst[order(rowSums(assay(vst)),decreasing=T),]

# output file
pdf(paste(RHB,"time_graphs_v2.pdf",sep="_"))

# plotting function
plotOTUs(assay(vst),vst@colData,facet=formula(Year + Country ~ OTU),line="smooth",design="time",colour="Treatment",plotsPerPage=6)

dev.off()

#===============================================================================
#       Yeast X2 sample
#===============================================================================
		    
# Only two samples, so can't do much more than descriptive statistics + graphs
dds2 <- dds[,dds$Treatment=="Yeast X 2"]
dds2 <- dds2[rowSums(counts(dds2,normalize=T))>5,]
output1 <- data.table(inner_join(data.table(OTU=rownames(dds2),W1=counts(dds2,normalize=T)[,1],W4=counts(dds2,normalize=T)[,2]),data.table(OTU=rownames(taxData),taxData)))
write.table(output1,paste(RHB,"Yx2_OTUs.txt",sep="_"),sep="\t",row.names=F,quote=F)

vst<-varianceStabilizingTransformation(dds)
vst <- vst[row.names(dds2),((dds$Time.point=="1 week")|(dds$Time.point=="4 week"))&dds$Country=="Germany"]
vst <- vst[order(rowSums(assay(vst)),decreasing=T),]

plotOTUs(assay(vst),X,facet=formula(~OTU),line="straight",design="time",colour="Treatment",plotsPerPage=4)
dev.off()

		     
