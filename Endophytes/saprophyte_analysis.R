#===============================================================================
#       Load libraries 
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(dplyr)
library(plyr)
library(Biostrings)
library(locfit)
library(ggplot2)
library(viridis)
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

# save data into a list, then ubiom_16S$countData to access countData and etc.
ubiom_16S <- list(
	countData=countData,
	colData=colData[colData$Loci!="18S",],
	taxData=taxData
)

# remove 18S level from colData
ubiom_16S$colData$Loci <- droplevels(ubiom_16S$colData$Loci)

# or all in one
ubiom_ITS <- list(
	countData=read.table("FUN.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=colData[colData$Loci!="16S",],
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=T,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65)
) 
ubiom_ITS$colData$Loci <- droplevels(ubiom_ITS$colData$Loci)

#===============================================================================
#       Create DEseq objects 
#===============================================================================

# "attach" the 16S verion of countData, colData and TaxData, three versions supplied - the third allows you to put the objects in their own namespace
# attach(ubiom_16S)
# apply(seq_along(ubiom_16S),function(i) assign(names(ubiom_16S)[i], ubiom_16S[[i]]))
 invisible(mapply(assign, names(ubiom_ITS), ubiom_ITS, MoreArgs=list(envir = globalenv())))


# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - using geoMeans function due to DES error
# if every gene contains at least one zero, cannot compute log geometric means
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)

# remove batch effect with limma removeBatchEffect method (not worth doing here)
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

# as above
d <- t(data.frame(t(pc.res)*mypca$percentVar))
# plot the PCA
pdf("saprophyte.its.pdf")
plotOrd(df,dds@colData,design="Country",shape="Year")
plotOrd(df,dds@colData,shape="Treatment",design="Time.point")
plotOrd(d,dds@colData,design="Country",shape="Year")
plotOrd(d,dds@colData,shape="Treatment",design="Time.point")
dev.off()

# saprophyte sample X205_S37 doesn't cluster as well as the other samples

#===============================================================================
#       Differential analysis of bean data
#===============================================================================

# get rid of the yeast2 samples (or the model will be unbalanced and "stuff")
dds <- dds[,(colnames(dds)!="X1_S88")&(colnames(dds)!="X4_S96")]
dds$Treatment <- droplevels(dds$Treatment)

# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

# the full model 
full_design <- ~Year + Country + Time.point + Treatment  + Treatment:Time.point

# the reduced model
reduced_design <- ~Year + Country + Treatment  + Time.point

design(dds) <- full_design

dds <-DESeq(dds, betaPrior=FALSE, test="LRT",full=full_design,reduced=reduced_design,parallel=T)

# calculate fit
dds <- DESeq(dds,parallel=T)

# Calculate OTUs which respond differently over time (time points) due to the treatment
res <- results(dds,parallel=T)

# 
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]


# contrasts to test
# condition effect
contrast=c("condition","Urea", "Control" )
contrast=c("condition","Yeast", "Control" )
contrast=c("condition","Urea", "Yeast" )

# calculate results
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)

# merge results with taxonomy table
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

write.table(res.merge,"saprophyte.its.urea_control.txt",quote=F,sep="\t",na="",row.names=F)
write.table(res.merge,"saprophyte.its.yeast_control.txt",quote=F,sep="\t",na="",row.names=F)
write.table(res.merge,"saprophyte.its.urea_yeast.txt",quote=F,sep="\t",na="",row.names=F)


# time effect
contrast=rbind( 
c("date","2_week", "1_week"),
c("date","4_week", "1_week"),
c("date","8_week", "1_week"),
c("date","16_week", "1_week"),
c("date","4_week", "2_week"),
c("date","8_week", "2_week"),
c("date","16_week", "2_week"),
c("date","8_week", "4_week"),
c("date","16_week", "4_week"),
c("date","16_week", "8_week"))

res <- apply(contrast,1,function(x) results(dds,alpha=alpha,parallel=T,contrast=x))
#sapply(res,summary)

res.merge <- lapply(res,function(o) data.table(inner_join(data.table(OTU=rownames(o),as.data.frame(o)),data.table(OTU=rownames(taxData),taxData))))

sapply(seq(1,10),function(i) write.table(res.merge[i],paste("Saprophyte.ITS",contrast[i,2],contrast[i,3],"txt",sep="."),quote=F,sep="\t",na="",row.names=F))


# interaction term (resultsNames(dds))
# contrast=list( "beanFRENCH.conditionHEALTHY","beanFrench.conditionHEALTHY") # effect of condition on french beans 
contrast=list( "beanRUNNER.conditionHEALTHY","beanRUNNER.conditionSICK") # effect of condition on runner beans

#===============================================================================
#       Graph analysis
#===============================================================================

dds2 <-dds[rowMeans(counts(dds,normalize=T))>=1,]
design(dds2) <- ~farm + date + condition

rld <- rlog(dds2,blind=F)
rld <- rld[order(rowSums(assay(rld)),decreasing=T),]
rld$colour <- paste(rld$farm,rld$condition,sep="_")
rld$time <- as.numeric(sub("_week","",rld$date))


rld2 <- log2(counts(dds,normalize=T)+1)
rld2 <- rld2[order(rowSums(rld2),decreasing=T),]

d <- data.frame(t(assay(rld[1,])),rld@colData)

g <- ggplot(data=d,aes_string(y=colnames(d)[1], x="time",colour="condition"))
g <- g + theme_classic_thin()
g <- g + scale_colour_viridis(discrete=TRUE)
g <- g + geom_point(size=2)
g <- g + facet_grid(.~ farm)
g + stat_smooth(method=locfit, formula=y~lp(x),se=F) 

d <- data.frame(t(assay(rld)),rld@colData)
d <- melt(d,id.vars = colnames(d)[(ncol(d)-10):ncol(d)],variable.name = "OTU", value.name = "rlog_counts")

d$rlog_counts <- d$rlog_counts+abs(min(d$rlog_counts))
ymax <- max(d$rlog_counts)


pdf("mega_wide.pdf",height=12,width=20)
noPlots <- 19
allVars <- unique(d$OTU)
noVars <- length(allVars)

plotSequence <- c(seq(0, noVars-1, by = noPlots), noVars)


sapply(seq(2,length(plotSequence)),function(i) {
	start <- plotSequence[i-1] + 1
	end <- plotSequence[i]
	tmp <- d[d$OTU %in% allVars[start:end],]
	cat(unique(tmp$OTU), "\n")

	g <- ggplot(data=tmp,aes(y=rlog_counts, x=time,colour=condition),ylim=c(0,ymax))
	g <- g + theme_classic_thin(base_size = 16) %+replace% theme(panel.border=element_rect(colour="black",size=0.25,fill=NA),legend.position="bottom")
	g <- g + scale_colour_viridis(discrete=TRUE)
	g <- g + facet_grid(farm ~ OTU,scales="free_x")
	g <- g + geom_point(size=2)
	g <- g + stat_smooth(method=locfit, formula=y~lp(x),se=F)
	print(g)
})

dev.off()



g <- ggplot(data=d,aes(y=counts, x=time,colour=condition))
g <- g + theme_classic_thin() %+replace% theme(panel.border=element_rect(colour="black",size=0.25,fill=NA))
g <- g + scale_colour_viridis(discrete=TRUE)
g <- g + geom_point(size=2)
g + stat_smooth(method=locfit, formula=y~lp(x),se=F) + facet_grid(OTU ~ farm)
g + stat_smooth(method=locfit, formula=y~lp(x),se=F) + facet_grid(farm ~ OTU,scales="free_y")

#===============================================================================
#       Differential OTUs
#===============================================================================

# read OTU fasta
myOTUs <- readDNAStringSet("ITS.otus.fa")

writeXStringSet(myOTUs[as.matrix(res.merge[padj<=0.1,1]),"main_big.fa")

#===============================================================================
#       PCA plots of beans only
#===============================================================================

mypca <- des_to_pca(dds)
df <-t(data.frame(t(mypca$x)*mypca$percentVar))

plotOrd(df,dds@colData,design="condition",title="Condition")
plotOrd(df,dds@colData,design="farm",title="Farm")
plotOrd(df,dds@colData,design="variety", shape="bean",title="Bean Type")
plotOrd(df,dds@colData,design="field", shape="variety",title="Field")



