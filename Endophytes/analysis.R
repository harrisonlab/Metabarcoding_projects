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
ubiom_BAC <- list(
	countData=countData,
	colData=colData[colData$Loci!="18S",],
	taxData=taxData,
	RHB="BAC"
)

# remove 18S level from colData
ubiom_BAC$colData$Loci <- droplevels(ubiom_BAC$colData$Loci)

# or all in one
ubiom_FUN <- list(
	countData=read.table("FUN.zotus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=colData[colData$Loci!="16S",],
	taxData=phyloTaxaTidy(read.table("zFUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 
ubiom_FUN$colData$Loci <- droplevels(ubiom_FUN$colData$Loci)

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

# calculate size factors - use geoMeans function if every gene contains at least one zero
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)


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
pdf(paste(RHB,"saprophyte.pdf",sep="_"))
plotOrd(df,dds@colData,design="Year",shape="Country")
plotOrd(d,dds@colData,design="Year",shape="Country")
plotOrd(df,dds@colData,shape="Treatment",design="Time.point")
plotOrd(d,dds@colData,shape="Treatment",design="Time.point")
dev.off()

# Fungal saprophyte sample X205_S37 doesn't cluster as well as the other samples

#===============================================================================
#       Differential analysis of saprophyte data
#===============================================================================

# get rid of the yeast2 samples (not enough samples to do any statistics with them)
dds <- dds[,(colnames(dds)!="X1_S88")&(colnames(dds)!="X4_S96")]

# drop the Yeast2 level from the Treatment factor
dds$Treatment <- droplevels(dds$Treatment)

# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power 
# but, no point keeping OTUs with 0 count
dds<-dds[rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

### Full model design ####

# the full model 
full_design <- ~Year + Country  + Treatment + Time.point + Treatment:Time.point

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# main effect urea vs control
contrast <- c("Treatment","Urea","Control")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Urea_effect.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# main effect yeast vs control
contrast <- c("Treatment","Yeast","Control")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Yeast_effect.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# treatment effect at each time point
contrast <- list("TreatmentYeast.Time.point1.week","TreatmentControl.Time.point1.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Yeast_W1.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentUrea.Time.point1.week","TreatmentControl.Time.point1.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Urea_W1.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentYeast.Time.point2.week","TreatmentControl.Time.point2.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Yeast_W2.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentUrea.Time.point2.week","TreatmentControl.Time.point2.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Urea_W2.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentYeast.Time.point4.week","TreatmentControl.Time.point4.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Yeast_W4.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentUrea.Time.point4.week","TreatmentControl.Time.point4.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Urea_W4.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentYeast.Time.point8.week","TreatmentControl.Time.point8.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Yeast_W8.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentUrea.Time.point8.week","TreatmentControl.Time.point8.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Urea_W8.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentYeast.Time.point16.week","TreatmentControl.Time.point16.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Yeast_W16.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

contrast <- list("TreatmentUrea.Time.point16.week","TreatmentControl.Time.point16.week")
res <-  results(dds,alpha=alpha,parallel=T,contrast=contrast)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"Urea_W16.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# a function for reading in and merging some/all the above files - base on a regex variable to specify which files
qfun <- function(regex_path){
	qq <- lapply(list.files(".",regex_path,full.names=T,recursive=F),function(x) fread(x)) # read in the files
	names(qq) <- list.files(".",regex_path,full.names=F,recursive=F) # gets the name of each file
	qq <- lapply(qq,function(l) l[,c(-4,-5,-6)]) # drops "lfcSE", "stat" and "pvalue" columns
	qq <- Map(function(x, i) {
		colnames(x)[3]<-paste(i, colnames(x)[3],sep="_");
		colnames(x)[4]<-paste(i, colnames(x)[4],sep="_");
		return(x)
	},qq, sub("\\.txt","",names(qq)))
	m <- Reduce(function(...) merge(..., all = T), qq) # could use inner_join rather than merge - would save maybe a couple of milliseconds
	return(m)
}
	
write.table(qfun("FUN_Urea.*.txt$"),"FUN_UREA_ALL.txt",sep="\t",row.name=F,quote=F)
write.table(qfun("FUN_Yeast.*.txt$"),"FUN_YEAST_ALL.txt",sep="\t",row.name=F,quote=F)
write.table(qfun("BAC_Urea.*.txt$"),"BAC_UREA_ALL.txt",sep="\t",row.name=F,quote=F)
write.table(qfun("BAC_Yeast.*.txt$"),"BAC_YEAST_ALL.txt",sep="\t",row.name=F,quote=F)


## difference over time ##

# the reduced model (for calculating rsponse to time)
reduced_design <- ~Year + Country + Treatment  + Time.point

# calculate model, including both full and reduced designs
dds <-DESeq(dds, betaPrior=FALSE, test="LRT",full=full_design,reduced=reduced_design,parallel=T)

# calculate OTUs which respond differently over time (time points) due to the treatment (this is the full vs reduced model)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))

# write table to file
write.table(res.merge, paste(RHB,"time_effect.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# this looks interesting, extract the results for each treatment?...
res <- results(dds,alpha=alpha,parallel=T,name= "Treatment_Urea_vs_Control",test="Wald")
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"time_effect_urea.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

res <- results(dds,alpha=alpha,parallel=T,name= "Treatment_Yeast_vs_Control",test="Wald")
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"time_effect_yeat.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)


write.table(qfun("FUN_time.*.txt$"),"FUN_TIME_ALL.txt",sep="\t",row.name=F,quote=F)
write.table(qfun("BAC_time.*.txt$"),"BAC_TIME_ALL.txt",sep="\t",row.name=F,quote=F)

### simplfied models ###

# 2017 data only

dds2 <- dds[,dds$Year==2017]
dds2$Treatment <- droplevels(dds2$Treatment)

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

#### test plot for first OTU
#d <- data.frame(t(assay(rld[1,])),rld@colData)
#d$time <- as.integer(sub(" week","",rld$Time.point))
#g <- ggplot(data=d,aes_string(y=colnames(d)[1], x="time",colour="Treatment"))
#g <- g + theme_classic_thin() %+replace% theme(panel.border=element_rect(colour="black",size=0.25,fill=NA),legend.position="bottom")
#g <- g + scale_colour_viridis(discrete=TRUE)
#g <- g + geom_point(size=2)
#g <- g + facet_grid(.~ Year + Country)
#g <- g + stat_smooth(method=locfit, formula=y~lp(x),se=F) 
#g
### end test plot

# extract rld values and combine with colData
d <- data.frame(t(assay(rld)),rld@colData)

# add an integer time column
d$time <- as.integer(sub(" week","",rld$Time.point))

# convert data to "long" format (single column will contain all the OTUs rather than one column per OTU) 
d <- melt(d,id.vars = colnames(d)[(ncol(d)-6):ncol(d)],variable.name = "OTU", value.name = "rlog_counts")

# zero bound rlog values (no negative values on y-axis here)
d$rlog_counts <- d$rlog_counts+abs(min(d$rlog_counts))

# set the maximum extent of the y-axis (all graphs will be on same scale)
ymax <- max(d$rlog_counts)

# number of plots per page
noPlots <- 5

# get unique OTUs
allVars <- unique(d$OTU)

# number of unique OTUs same as: noVars <- length(unique(d$OTU)) 
noVars <- length(allVars)

# plot holder
plotSequence <- c(seq(0, noVars-1, by = noPlots), noVars)

# output file
pdf(paste(RHB,"time_graphs.pdf",sep="_"))

# plotting function
sapply(seq(2,length(plotSequence)),function(i) {
	start <- plotSequence[i-1] + 1
	end <- plotSequence[i]
	tmp <- d[d$OTU %in% allVars[start:end],]
	cat(unique(tmp$OTU), "\n")
	g <- ggplot(data=tmp,aes(y=rlog_counts, x=time,colour=Treatment),ylim=c(0,ymax))
	g <- g + theme_classic_thin(base_size = 16) %+replace% theme(panel.border=element_rect(colour="black",size=0.25,fill=NA),legend.position="bottom")
	g <- g + scale_colour_viridis(discrete=TRUE)
	g <- g + facet_grid(Year + Country ~ OTU,scales="free_x")
	g <- g + geom_point(size=2)
	g <- g + stat_smooth(method=locfit, formula=y~lp(x),se=F)
	print(g)
})

dev.off()
