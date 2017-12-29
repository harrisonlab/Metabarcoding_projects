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
library(Biostrings)

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

# save data into a list
ubiom_BAC <- list(
	countData=countData,
	colData=colData,
	taxData=taxData,
	RHB="BAC"
)

# Fungi all in one call
ubiom_FUN <- list(
	countData=read.table("FUN.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="FUN"
) 

# Oomycetes
ubiom_OO <- list(
	countData=read.table("OO.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData2",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("OO.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="OO"
) 
rownames(ubiom_OO$colData) <- paste0("X",gsub("_","\\.",ubiom_OO$colData$name),"_",sub("D.*","",rownames(ubiom_OO$colData)))

# Nematodes 
ubiom_NEM <- list(
	countData=read.table("NEM.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("colData2",header=T,sep="\t",row.names=1),
	taxData=phyloTaxaTidy(read.table("NEM.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
	RHB="NEM"
) 
rownames(ubiom_NEM$colData) <- paste0("X",gsub("_","\\.",ubiom_NEM$colData$name),"_",sub("D.*","",rownames(ubiom_NEM$colData)))

#===============================================================================
#       Combine species 
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species) - no oo or nem in this data set

# oomycetes

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("FUN.taxa")
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData

# Nematodes

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
#rownames(colData) <- sub("^XG","G",rownames(colData))
colData <- colData[names(countData),]

# remove low count samples and control samples (not needed here)
filter <- (colSums(countData)>=1000)&colData$condition!="C"
colData <- droplevels(colData[filter,])
countData <- countData[,filter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - use geoMeans function if
# every gene contains at least one zero, as cannot compute log geometric means
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# sizeFactors(dds) <-geoMeans(dds)
# library(edgeR) # I think anyway
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

# Correction from aboslute quantification
sizeFactors(dds) <- sizeFactors(dds)/colData$funq
sizeFactors(dds) <- sizeFactors(dds)/colData$bacq

#===============================================================================
#       Filter data 
#============================================================================

### filters to remove OTUs which are unlikely part of the correct kingdom (SAR and 18S Eukaryote)
# pythium
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]
# nematode
myfilter <- row.names(taxData[as.number(taxData$c_conf)>0.9 & as.number(taxData$o_conf)>0.9,])
dds <- dds[rownames(dds)%in%myfilter,]

### read accumulation filter
# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

#### Select filter ####
# Apply seperately for appropriate data set depending on cut-off chosen from graph
myfilter <- dtt$OTU[1:500] #FUN
myfilter <- dtt$OTU[1:50] # OO
myfilter <- dtt$OTU[1:20] # NEM
myfilter <- dtt$OTU[1:1000]  # BAC

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
colData$location<-as.number(colData$pair)

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
plotOrd(df,colData,design="condition",xlabel="PC1",ylabel="PC2")
plotOrd(df,colData,shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

### remove spatial information (this uses the factor "pair" not the numeric "location") and plot
pc.res <- resid(aov(mypca$x~colData$pair,colData))
df <- t(data.frame(t(pc.res*mypca$percentVar)))

pdf(paste(RHB,"PCA_deloc.pdf",sep="_"))
plotOrd(df,colData,shape="condition",design="location",continuous=T,xlabel="PC1",ylabel="PC2")
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

# contrast (not actually necessary in this case as this would be the default result calculated by results(dds)
# contrast <- c("condition","S","H")
res <- results(dds,alpha=alpha,parallel=T)
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.1]$OTU],paste0(RHB,".sig.fa"))


#===============================================================================
#       Pie chart
#===============================================================================
taxData<-phyloTaxaTidy(taxData,0.9,level=2)

X <- sumTaxaAdvanced(list(as.data.frame(counts(dds,normalized=F)),taxData,colData(dds)),taxon="rank",proportional=T,cutoff=1)
Y <- sumTaxaAdvanced(list(as.data.frame(counts(dds,normalized=T)),taxData,colData(dds)),taxon="rank",proportional=T,cutoff=1)
X$type=1
Y$type=2
X <- rbind(X,Y)

X$Phylum <- droplevels(X$rank)
levels(X$Phylum) <- sub("\\(.*","",X$Phylum)
#levels(X$Phylum)[5] <- "Candidatus Saccharibacteria"

X <- X %>% group_by(type) %>% mutate(pos = 100-(cumsum(all)- all/2))

g <- ggplot(X,aes(x="",y=all,fill=Phylum,label = sprintf("%0.1f", round(all, digits = 1))))
g <- g + geom_bar(width = 1,stat="identity")
g <- g + scale_fill_viridis(discrete=TRUE,direction=-1,begin=0.2)
g <- g + geom_text(data=X[X$all>5,],aes(y = pos),size=3)
#g <- g +geom_text(aes(label = sprintf("%0.2f", round(all, digits = 2))), position = position_stack(vjust = 0.5,),size=2)
g <- g + geom_label_repel(data=X[X$all<=5,],aes(y = pos), size=3, show.legend = F, nudge_x=0.75,nudge_y=0.75,segment.size=0.1)
g <- g + coord_polar(theta="y",direction=-1)
g <- g +facet_grid(facets=. ~ type,switch="x")
g <- g + theme_blank() %+replace% theme(	
#  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  axis.text  = element_blank(),
  axis.ticks = element_blank(),
  strip.background = element_blank(),
  #strip.text.x = element_blank()
)
ggsave(paste0(RHB,"_pie.pdf"),g+xlab(NULL)+ylab(NULL))


#===============================================================================

#       Alpha diversity analysis

#===============================================================================



# Recreate dds object and don't filter for low counts before running Alpha diversity



# add some extra columns to coldata to describe the data (easier than releveling condition and changing row names)

colData$Sample<-c(paste0(rep("B",3),seq(1,3)),paste0(rep("E",3),seq(1,3)),paste0(rep("A",3),seq(1,3)),paste0(rep("P",3),seq(1,3)))

colData$type <- as.factor(c(rep("Between",3),rep("Edge",3),rep("Apple",3),rep("Pear",3)))



# plot alpha diversity - plot_alpha will convert normalised abundances to integer values (limits for bac only)

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="type",colour="Sample"))



### permutation based anova on diversity index ranks ###



# get the diversity index data

all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="condition",colour="block",returnData=T)



# add column names as row to metadata (or use tribble)

colData$Samples <- rownames(colData)



# join diversity indices and metadata

all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData))



# perform anova for each index

sink(paste(RHB,"ALPHA_stats.txt",sep="_"))

setkey(all_alpha_ord,S.chao1)

summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~condition+Error(treatment),all_alpha_ord))

setkey(all_alpha_ord,shannon)

summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~condition+Error(treatment),all_alpha_ord))

setkey(all_alpha_ord,simpson)

summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~condition+Error(treatment),all_alpha_ord))

setkey(all_alpha_ord,S.ACE)

summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~condition+Error(treatment),all_alpha_ord))

sink()

#===============================================================================

#       ANOVA

#===============================================================================



# Want to see if genotype can describe more of the variance in the rhizosphere data compared to the bulk soil data.

# This would suggest that the rhizospere are recruting specific OTUs



qfun <- function (mypca,dds,names=c("area","genotype","interaction","residual"),model) {

	X <- t(apply(mypca$x,2,function(x){t(summary(aov(as.formula(paste0("x",model)),colData(dds)))[[1]][2])}))

	colnames(X) <- names

	x<-t(apply(X,1,prop.table))

	perVar <- x * mypca$percentVar

	print(colSums(perVar))

	print(colSums(perVar)/sum(colSums(perVar))*100)

	return(as.data.table(perVar))

}



#pc.res <- resid(aov(mypca$x~run,colData(dds)))



#ss_all_res <- qfun(list(x=pc.res,percentVar=mypca$percentVar),dds)



ss_all <- qfun(mypca,dds,model="~area + genotype + area * genotype",names=c("area","genotype","interaction","residual"))

ss_all_run <- qfun(mypca,dds,model="~run+area + genotype + area * genotype",names=c("run","area","genotype","interaction","residual"))

lapply(seq(1:4),function(x) {X<-summary(aov(mypca$x[,x]~area + genotype + area * genotype,colData(dds)))[[1]][[2]];X/sum(X)*100})

lapply(seq(1:4),function(x) {X<-summary(aov(mypca$x[,x]~run+area + genotype + area * genotype,colData(dds)))[[1]][[2]];X/sum(X)*100})



dds_rhiz  <- dds[,dds$area=="Rhizosphere"]

dds_stool <- dds[,dds$area!="Rhizosphere"]



mypca_rhiz <- des_to_pca(dds_rhiz)

mypca_stool <- des_to_pca(dds_stool)



# sum of squares

ss_pervar_rhiz  <- qfun(mypca_rhiz,dds_rhiz,names=c("run","genotype","residual"),model="~run+genotype")

ss_pervar_rhiz  <- qfun(mypca_rhiz,dds_rhiz,names=c("genotype","residual"),model="~genotype")

ss_pervar_stool <- qfun(mypca_stool,dds_stool,names=c("run","genotype","residual"),model="~run+genotype")

ss_pervar_stool <- qfun(mypca_stool,dds_stool,names=c("genotype","residual"),model="~genotype")



lapply(seq(1:4),function(x) {X<-summary(aov(mypca_rhiz$x[,x]~run+genotype,colData(dds_rhiz)))[[1]][[2]];X/sum(X)*100})

lapply(seq(1:4),function(x) {X<-summary(aov(mypca_rhiz$x[,x]~genotype,colData(dds_rhiz)))[[1]][[2]];X/sum(X)*100})



lapply(seq(1:4),function(x) {X<-summary(aov(mypca_stool$x[,x]~run+genotype,colData(dds_stool)))[[1]][[2]];X/sum(X)*100})

lapply(seq(1:4),function(x) {X<-summary(aov(mypca_stool$x[,x]~genotype,colData(dds_stool)))[[1]][[2]];X/sum(X)*100})



# alternative - combine both P26 sand and clay into single genotype

levels(dds$genotype)[2:3] <- "M26"

mypca <- des_to_pca(dds)

# then re-run this section


