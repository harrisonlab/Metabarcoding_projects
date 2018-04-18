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
library(outliers)

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

ubiom_FUN$taxData <- ubiom_FUN$taxData[rownames(ubiom_FUN$countData),]
ubiom_BAC$taxData <- ubiom_BAC$taxData[rownames(ubiom_BAC$countData),]

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
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC,calcFactors=geoMeans)

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Pool Data/subsample
#===============================================================================

# Get rid of control samples
# could be interesting to see if the control samples are less different over time than the treated samples
dds <- dds[,gsub("(^[A-Z][0-9]*)([A-Z])(.*)","\\2",rownames(colData(dds)))!="C"]

# There are only 3 (out of 900) missing samples - subsampling is a bit too extreme
# As each sample point has three biological replicates will take the mean of other samples to represent missing samples

# get number of samples per tree
sample_numbers <- table(sub("[A-Z]$","",rownames(colData(dds))))

# collapse (sum) samples
dds <- collapseReplicates2(dds,groupby=sub("[A-Z]$","",rownames(colData(dds))),simple=T)

# set the dds sizefactor to the number of samples
dds$sizeFactor <- as.vector(sample_numbers/3)

# recreate countData and colData
countData<- round(counts(dds,normalize=T),0)
colData <- as.data.frame(colData(dds))

# new dds object with the corrected data set
dds <- DESeqDataSetFromMatrix(countData,colData,~1)

# calculate size factors - using geoMeans function (works better with this data set)
max(geoMeans(dds))/min(geoMeans(dds))
max(sizeFactors(estimateSizeFactors(dds)))/min(sizeFactors(estimateSizeFactors(dds)))
sizeFactors(dds) <-geoMeans(dds)
# sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

dds$site 			 <- substr(colnames(dds),1,1)
dds$loc_factor <- as.factor(dds$meters)
dds$time  		 <- as.factor(dds$time)
dds$block			 <- as.factor(dds$block)

list_dds <-list(all     = dds,
		cider   = dds[,dds$site=="H"],
		dessert = dds[,dds$site=="G"],
		c_tree  = dds[,dds$site=="H"&dds$condition=="Y"],
		c_grass = dds[,dds$site=="H"&dds$condition=="N"],
		d_tree  = dds[,dds$site=="G"&dds$condition=="Y"],
		d_grass = dds[,dds$site=="G"&dds$condition=="N"])

#===============================================================================
#       Alpha diversity analysis - RUN BEFORE FILTERING OUT ANY LOW COUNT OTUS
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

ggsave(paste(RHB,"Alpha_all.pdf",sep="_"),plot_alpha(counts(list_dds$all,normalize=T),colData(list_dds$all),colour="site",design="time",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))
ggsave(paste(RHB,"Alpha_cider.pdf",sep="_"),plot_alpha(counts(list_dds$cider,normalize=T),colData(list_dds$cider),colour="condition",design="time",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))
ggsave(paste(RHB,"Alpha_dessert.pdf",sep="_"),plot_alpha(counts(list_dds$dessert,normalize=T),colData(list_dds$dessert),colour="condition",design="time",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))
ggsave(paste(RHB,"Alpha_c_tree.pdf",sep="_"),plot_alpha(counts(list_dds$c_tree,normalize=T),colData(list_dds$c_tree),design="genotype_name",colour="time",discrete=T,measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))
ggsave(paste(RHB,"Alpha_c_grass.pdf",sep="_"),plot_alpha(counts(list_dds$c_grass,normalize=T),colData(list_dds$c_grass),design="genotype_name",colour="time",discrete=T,measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))
ggsave(paste(RHB,"Alpha_d_tree.pdf",sep="_"),plot_alpha(counts(list_dds$d_tree,normalize=T),colData(list_dds$d_tree),design="genotype_name",colour="time",discrete=T,measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))
ggsave(paste(RHB,"Alpha_d_grass.pdf",sep="_"),plot_alpha(counts(list_dds$d_grass,normalize=T),colData(list_dds$d_grass),design="genotype_name",colour="time",discrete=T,measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))

### permutation based anova on diversity index ranks ###

# get alpha diversity indices
list_alpha_ord <- lapply(list_dds,function(dds) plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",returnData=T))

# join diversity indices and metadata
list_alpha_ord <- lapply(seq_along(list_dds),function(i) {
	as.data.table(left_join(list_alpha_ord[[i]],as.data.frame(colData(list_dds[[i]]))%>% mutate(Samples = rownames(colData(list_dds[[i]])))))
})


# perform anova for each index  - need to work out what we need to test for

sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
  all_alpha_ord <- list_alpha_ord[[1]]
	setkey(all_alpha_ord,S.chao1)
	print("Chao1")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~site+loc_factor+condition*time*genotype_name,all_alpha_ord))
	setkey(all_alpha_ord,shannon)
	print("Shannon")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~site+loc_factor+condition*time*genotype_name,all_alpha_ord))
	setkey(all_alpha_ord,simpson)
	print("simpson")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~site+loc_factor+condition*time*genotype_name,all_alpha_ord))
sink()

summary(aov(list_pca[[1]]$x[,x]~site+loc_factor+condition+time*genotype_name,colData(list_dds[[1]])))


sink(paste(RHB,"ALPHA_stats_sep.txt",sep="_"))
	lapply(list_alpha_ord[4:7],function(all_alpha_ord) {
	  setkey(all_alpha_ord,S.chao1)
	  print("Chao1")
	  print(summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~loc_factor+time*genotype_name,all_alpha_ord)))
  	setkey(all_alpha_ord,shannon)
  	print("Shannon")
	  print(summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~loc_factor+time*genotype_name,all_alpha_ord)))
	  setkey(all_alpha_ord,simpson)
	  print("simpson")
	  print(summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~loc_factor+time*genotype_name,all_alpha_ord)))
	})
sink()


#===============================================================================
#       Filter data
#===============================================================================

# filter OTUs
otufilter <- lapply(list_dds,function(dds) rowSums(counts(dds, normalize=T))>4)

# filter count data
list_dds <- lapply(seq_along(list_dds),function(i) list_dds[[i]][otufilter[[i]],])

# filter taxonomy data
list_taxData <- lapply(seq_along(list_dds),function(i) taxData[rownames(list_dds[[i]]),])

#===============================================================================
#       Microbial Populations
#===============================================================================

### phylum level population frequency ###

#Dessert orchard
sink("Phylum_Frequencies.txt")
 cat("# Dessert Orchard frequenciea at phylum rank\n")
 cat("# Overall\n")
 sumTaxa(list(as.data.frame(counts(list_dds$dessert,normalize=T)),taxData,list_dds$dessert@colData),conf=0.8,design="all",proportional=T)
 cat("# By time point\n")
 sumTaxa(list(as.data.frame(counts(list_dds$dessert,normalize=T)),taxData,list_dds$dessert@colData),conf=0.8,design="time",proportional=T)
 cat("# By genotype\n")
 sumTaxa(list(as.data.frame(counts(list_dds$dessert,normalize=T)),taxData,list_dds$dessert@colData),conf=0.8,design="genotype_name",proportional=T)
 cat("# By time point and condition\n")
 sumTaxa(list(as.data.frame(counts(list_dds$dessert,normalize=T)),taxData,list_dds$dessert@colData),conf=0.8,design=c("condition","time"),proportional=T)
 cat("# By time point,condition and genotype\n")
 sumTaxa(list(as.data.frame(counts(list_dds$dessert,normalize=T)),taxData,list_dds$dessert@colData),conf=0.8,design=c("genotype_name","condition","time"),proportional=T)
 cat("\n\n# Cider Orchard frequenciea at phylum rank\n")
 cat("# Overall\n")
 sumTaxa(list(as.data.frame(counts(list_dds$cider,normalize=T)),taxData,list_dds$cider@colData),conf=0.8,design="all",proportional=T)
 cat("# By time point\n")
 sumTaxa(list(as.data.frame(counts(list_dds$cider,normalize=T)),taxData,list_dds$cider@colData),conf=0.8,design="time",proportional=T)
 cat("# By genotype\n")
 sumTaxa(list(as.data.frame(counts(list_dds$cider,normalize=T)),taxData,list_dds$cider@colData),conf=0.8,design="genotype_name",proportional=T)
 cat("# By time point and condition\n")
 sumTaxa(list(as.data.frame(counts(list_dds$cider,normalize=T)),taxData,list_dds$cider@colData),conf=0.8,design=c("condition","time"),proportional=T)
 cat("# By time point,condition and genotype\n")
 sumTaxa(list(as.data.frame(counts(list_dds$cider,normalize=T)),taxData,list_dds$cider@colData),conf=0.8,design=c("genotype_name","condition","time"),proportional=T)
sink()


## Run seperatly for each orchard ##
dds <- list_dds$dessert
ORCH<-"Dessert"

#dds <- list_dds$cider
#ORCH<-"Cider"

md <- melt(sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,dds@colData),conf=0.8,design=c("genotype_name","condition","time"),proportional=T))
md$Genotype  <-  unlist(strsplit(as.character(md[,2])," : "))[c(T,F,F)]
md$Condition <-  unlist(strsplit(as.character(md[,2])," : "))[c(F,T,F)]
md$TimePoint <-  unlist(strsplit(as.character(md[,2])," : "))[c(F,F,T)]

md$Condition[md$Condition=="N"] <- "Grass Alley"
md$Condition[md$Condition=="Y"] <- "Tree Station"
md$Condition[md$Condition=="C"] <- "Control"

md$Genotype  <- as.factor(md$Genotype)
md$Condition <- as.factor(md$Condition)
md$TimePoint <- as.integer(md$TimePoint)

md$phylum <- factor(md$phylum,aggregate(md$value,by=list(md$phylum),mean)[order(aggregate(md$value,by=list(md$phylum),mean)[,2],decreasing=T),][,1])

# plot
g <- ggplot(md,aes(x=TimePoint,y=value,colour=Condition,phylum=phylum))
g <- g + geom_smooth() + scale_x_continuous(breaks=c(0,1,2))
g <- g + scale_colour_manual(values = c("Grass Alley" = "black", "Tree Station" = "orange"))
g <- g + facet_wrap(~phylum,scales="free_y")
g <- g +  theme_classic_thin() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave(paste(RHB,ORCH,"Proportional graphs.pdf",sep="_"),g)

md <- melt(sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,dds@colData),conf=0.8,design=c("genotype_name","condition","time"),proportional=F,meanDiff=T))
md$Genotype  <-  unlist(strsplit(as.character(md[,2]),"\\.+"))[c(T,F,F,F,F,F)]
md$Condition <-  unlist(strsplit(as.character(md[,2]),"\\.+"))[c(F,T,F,F,F,F)]
md$TimePoint <-  unlist(strsplit(as.character(md[,2]),"\\.+"))[c(F,F,T,F,F,F)]
md$Genotype_2  <-  unlist(strsplit(as.character(md[,2]),"\\.+"))[c(F,F,F,T,F,F)]
md$Condition_2 <-  unlist(strsplit(as.character(md[,2]),"\\.+"))[c(F,F,F,F,T,F)]
md$TimePoint_2 <-  unlist(strsplit(as.character(md[,2]),"\\.+"))[c(F,F,F,F,F,T)]
md <- md[(md$Genotype==md$Genotype_2)&(md$TimePoint==md$TimePoint_2),]

md$tp_size <- 0
md$tp_size[md$TimePoint==0] <- sum(md$value[md$TimePoint==0])
md$tp_size[md$TimePoint==1] <- sum(md$value[md$TimePoint==1])
md$tp_size[md$TimePoint==2] <- sum(md$value[md$TimePoint==2])
md$prop <- (md$value/md$tp_size)*100


md$Condition[md$Condition=="N"] <- "Grass Alley"
md$Condition[md$Condition=="Y"] <- "Tree Station"
md$Condition[md$Condition=="C"] <- "Control"

md$Genotype  <- as.factor(md$Genotype)
md$Condition <- as.factor(md$Condition)
md$TimePoint <- as.integer(md$TimePoint)

colnames(md)[c(3,11)] <- c("Absolute","Relative")

pdf(paste(RHB,ORCH,"Mean Difference abundance.pdf",sep="_"),width=8,height=7)
 g <- ggplot(md,aes(x=TimePoint,phylum=phylum))
 g <- g + geom_smooth(aes(y=Absolute))
 g <- g + scale_x_continuous(breaks=c(0,1,2))
 g <- g + facet_wrap(~phylum,scales="free_y")
 g <- g +  theme_classic_thin() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
 g;g %>% remove_geom("smooth",1,"y") + geom_smooth(aes(y=Relative))
dev.off()

pdf(paste(RHB,ORCH,"Genotype_mean_difference.pdf"),width=20,height=20)
 g <- g + facet_wrap(Genotype~phylum,scales="free_y")
 g;g %>% remove_geom("smooth",1,"y") + geom_smooth(aes(y=Relative))
dev.off()
#===============================================================================
#       Beta diversity analysis
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
list_pca <- lapply(list_dds,des_to_pca)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <- lapply(list_pca,function(mypca) t(data.frame(t(mypca$x)*mypca$percentVar)))

pc.res <- lapply(seq_along(list_pca),function(i) resid(aov(list_pca[[i]]$x~list_dds[[i]]$loc_factor)))
dd <- lapply(seq_along(list_pca),function(i) t(data.frame(t(pc.res[[i]])*list_pca[[i]]$percentVar)))

# plot the PCA - need to think about model
pdf(paste(RHB,"PCA.pdf",sep="_"))
 plotOrd(d[[1]],colData(list_dds[[1]]),design="site",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Both orchards")
 plotOrd(d[[2]],colData(list_dds[[2]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider orchard")
 plotOrd(d[[3]],colData(list_dds[[3]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert orchard")
 plotOrd(d[[4]],colData(list_dds[[4]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider tree station")
 plotOrd(d[[5]],colData(list_dds[[5]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider grass alley")
 plotOrd(d[[6]],colData(list_dds[[6]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert tree station")
 plotOrd(d[[7]],colData(list_dds[[7]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert grass alley")
dev.off()

pdf(paste(RHB,"PCA_LOCATION.pdf",sep="_"))
 plotOrd(dd[[1]],colData(list_dds[[1]]),design="site",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Both orchards")
 plotOrd(dd[[2]],colData(list_dds[[2]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider orchard")
 plotOrd(dd[[3]],colData(list_dds[[3]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert orchard")
 plotOrd(dd[[4]],colData(list_dds[[4]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider tree station")
 plotOrd(dd[[5]],colData(list_dds[[5]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider grass alley")
 plotOrd(dd[[6]],colData(list_dds[[6]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert tree station")
 plotOrd(dd[[7]],colData(list_dds[[7]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert grass alley")
dev.off()

# ANOVA
#sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
#	print("ANOVA")
#	lapply(seq(1:4),function(x) {
#		summary(aov(list_pca[[1]]$x[,x]~site+loc_factor+condition*time*genotype_name,colData(list_dds[[1]])))
#	})
#	print("PERMANOVA")
#	lapply(seq(1:4),function(x) {
#		summary(aovp(list_pca[[1]]$x[,x]~site+loc_factor+condition*time*genotype_name,colData(list_dds[[1]])))
#	})
#sink()

#colData(list_dds)$block <- as.factor(colData(list_dds)$block
#colData(list_dds)$block <- as.factor(colData(list_dds)$block

# overall anova scores (% sum of squares)

# get the model sum of squares for all pc scores
sum_squares_cider <- #cbind(
	t(apply(list_pca[[2]]$x,2,function(x) t(summary(aov(x~block+time*condition*genotype_name,colData(list_dds[[2]])))[[1]][2])))#,
	#t(apply(list_pca[[2]]$x,2,function(x) t(summary(aov(x~genotype_name,colData(list_dds[[2]])))[[1]][2])))[,1]
#)
#colnames(sum_squares_cider) <- c("location","time","condition","time:condition","time:genotype","condition:genotype","time:condition:genotype_name","residual","genotype")
colnames(sum_squares_cider) <- c("location","time","condition","genotype","time:condition","time:genotype","condition:genotype","time:condition:genotype","residual")

sum_squares_dessert <- #cbind(
	t(apply(list_pca[[3]]$x,2,function(x) t(summary(aov(x~block+time*condition*genotype_name,colData(list_dds[[3]])))[[1]][2])))#,
	#t(apply(list_pca[[3]]$x,2,function(x) t(summary(aov(x~genotype_name,colData(list_dds[[3]])))[[1]][2])))[,1]
#)
colnames(sum_squares_cider) <- c("location","time","condition","genotype","time:condition","time:genotype","condition:genotype","time:condition:genotype","residual")

# output sum of squares %
sink(paste(RHB,"PCA_scores_by_orchard(blocks).txt",sep="_"))
 perVar<-t(apply(sum_squares_cider,1,prop.table))* list_pca[[2]]$percentVar
 cat("# Cider\n")
 colSums(perVar)/sum(colSums(perVar))*100
 perVar<-t(apply(sum_squares_dessert,1,prop.table))* list_pca[[2]]$percentVar
 cat(" # Dessert\n")
 colSums(perVar)/sum(colSums(perVar))*100
sink()


# sum squares at each time point
sum_squares_cider_0 <- #cbind(
	t(apply(list_pca[[2]]$x[colData(list_dds[[2]])$time==0,],2,function(x)t(summary(aov(x~block+condition*genotype_name,colData(list_dds[[2]])[colData(list_dds[[2]])$time==0,]))[[1]][2])))#,
	#t(apply(list_pca[[2]]$x[colData(list_dds[[2]])$time==0,],2,function(x)t(summary(aov(x~genotype_name,colData(list_dds[[2]])[colData(list_dds[[2]])$time==0,]))[[1]][2])))[,1]
#)
sum_squares_cider_1 <- #cbind(
	t(apply(list_pca[[2]]$x[colData(list_dds[[2]])$time==1,],2,function(x)t(summary(aov(x~block+condition*genotype_name,colData(list_dds[[2]])[colData(list_dds[[2]])$time==1,]))[[1]][2])))#,
	#t(apply(list_pca[[2]]$x[colData(list_dds[[2]])$time==1,],2,function(x)t(summary(aov(x~genotype_name,colData(list_dds[[2]])[colData(list_dds[[2]])$time==1,]))[[1]][2])))[,1]
#)
sum_squares_cider_2 <- #cbind(
	t(apply(list_pca[[2]]$x[colData(list_dds[[2]])$time==2,],2,function(x)t(summary(aov(x~block+condition*genotype_name,colData(list_dds[[2]])[colData(list_dds[[2]])$time==2,]))[[1]][2])))#,
	#t(apply(list_pca[[2]]$x[colData(list_dds[[2]])$time==2,],2,function(x)t(summary(aov(x~genotype_name,colData(list_dds[[2]])[colData(list_dds[[2]])$time==2,]))[[1]][2])))[,1]
#)
sum_squares_dessert_0 <- #cbind(
	t(apply(list_pca[[3]]$x[colData(list_dds[[3]])$time==0,],2,function(x)t(summary(aov(x~block+condition*genotype_name,colData(list_dds[[3]])[colData(list_dds[[3]])$time==0,]))[[1]][2])))#,
	#t(apply(list_pca[[3]]$x[colData(list_dds[[3]])$time==0,],2,function(x)t(summary(aov(x~genotype_name,colData(list_dds[[3]])[colData(list_dds[[3]])$time==0,]))[[1]][2])))[,1]
#)
sum_squares_dessert_1 <- #cbind(
	t(apply(list_pca[[3]]$x[colData(list_dds[[3]])$time==1,],2,function(x)t(summary(aov(x~block+condition*genotype_name,colData(list_dds[[3]])[colData(list_dds[[3]])$time==1,]))[[1]][2])))#,
#	t(apply(list_pca[[3]]$x[colData(list_dds[[3]])$time==1,],2,function(x)t(summary(aov(x~genotype_name,colData(list_dds[[3]])[colData(list_dds[[3]])$time==1,]))[[1]][2])))[,1]
#)
sum_squares_dessert_2 <- #cbind(
	t(apply(list_pca[[3]]$x[colData(list_dds[[3]])$time==2,],2,function(x)t(summary(aov(x~block+condition*genotype_name,colData(list_dds[[3]])[colData(list_dds[[3]])$time==2,]))[[1]][2])))#,
	#t(apply(list_pca[[3]]$x[colData(list_dds[[3]])$time==2,],2,function(x)t(summary(aov(x~genotype_name,colData(list_dds[[3]])[colData(list_dds[[3]])$time==2,]))[[1]][2])))[,1]
#)
# set column names
colnames(sum_squares_cider_0)   <- colnames(sum_squares_cider_1)   <- colnames(sum_squares_cider_2)   <- c("location","condition","genotype","condition:genotype","residual")
colnames(sum_squares_dessert_0) <- colnames(sum_squares_dessert_1) <- colnames(sum_squares_dessert_2) <- c("location","condition","genotype","condition:genotype","residual")

# output sum of squares %
sink(paste(RHB,"PCA_scores_by_time(blocks).txt",sep="_"))
 cat("# Cider\n")
 cat("# Time 0\n")
 perVar<-t(apply(sum_squares_cider_0,1,prop.table))* list_pca[[2]]$percentVar[colData(list_dds[[2]])$time==0]
 colSums(perVar)/sum(colSums(perVar))*100
 cat("# Time 1\n")
 perVar<-t(apply(sum_squares_cider_1,1,prop.table))* list_pca[[2]]$percentVar[colData(list_dds[[2]])$time==1]
 colSums(perVar)/sum(colSums(perVar))*100
 cat("# Time 2\n")
 perVar<-t(apply(sum_squares_cider_2,1,prop.table))* list_pca[[2]]$percentVar[colData(list_dds[[2]])$time==2]
 colSums(perVar)/sum(colSums(perVar))*100
 cat("\n# Dessert\n")
 cat("# Time 0\n")
 perVar<-t(apply(sum_squares_dessert_0,1,prop.table))* list_pca[[2]]$percentVar[colData(list_dds[[3]])$time==0]
 colSums(perVar)/sum(colSums(perVar))*100
 cat("# Time 1\n")
 perVar<-t(apply(sum_squares_dessert_1,1,prop.table))* list_pca[[2]]$percentVar[colData(list_dds[[3]])$time==1]
 colSums(perVar)/sum(colSums(perVar))*100
 cat("# Time 2\n")
 perVar<-t(apply(sum_squares_dessert_2,1,prop.table))* list_pca[[2]]$percentVar[colData(list_dds[[3]])$time==2]
 colSums(perVar)/sum(colSums(perVar))*100
sink()

sink(paste(RHB,"PCA_ANOVA_by_orchard.txt",sep="_"))
	lapply(seq(2,3),function(i) {
	  print(names(list_dds)[i])
	  lapply(seq(1,4),function(x) {
		  print(summary(aov(list_pca[[i]]$x[,x]~loc_factor+time*condition*genotype_name,colData(list_dds[[i]]))))
		
ord_rda <- list(ord_rda1[[1]],ord_rda1[[2]][[1]],ord_rda1[[2]][[2]],ord_rda1[[3]][[1]],ord_rda1[[3]][[2]],ord_rda1[[3]][[3]],ord_rda1[[3]][[4]])

theme_set(theme_classic_thin())
pdf(paste(RHB,"RDA.pdf",sep="_"))
lapply(seq_along(ord_rda),function(i) plot_ordination(list_phylo[[i]], ord_rda[[i]], "samples", color="genotype_name",shape="time",title=names(list_dds)[i])+ geom_point(size=3.5,alpha=0.75))
dev.off()

# permutation anova of rda1 and rda 2
aov_rda <- lapply(ord_rda,function(o) anova.cca(o,permuations=999,by="terms",model="reduced",parallel=12))

sink(paste(RHB,"RDA_permutation_anova",sep="_"))
	aov_rda
sink()

#===============================================================================
#       RDA plots
#===============================================================================

# WORK OUT WHAT I WANT TO PLOT FIRST BEFORE WRITING THE CODE!!!!!




# agregate counts at phylum level
combinedTaxa <- lapply(list_taxData,function(taxData) combineTaxa2(taxData[,-8],rank="phylum",confidence=0.8,returnFull=T))
phylumData <- lapply(seq_along(combinedTaxa),function(i) combCounts(combinedTaxa[[i]],counts(list_dds[[i]],normalize=T)))
phylumTaxa <- lapply(seq_along(combinedTaxa),function(i) combTaxa(combinedTaxa[[i]],list_taxData[[i]][,-8]))
phylumData <- lapply(seq_along(phylumData), function(i) aggregate(phylumData[[i]],by=list(taxaConfVec(phylumTaxa[[i]],conf=0.8,level=2)),sum))
phylumData <- lapply(phylumData,function(phylumData) {rownames(phylumData) <- phylumData[,1];phylumData <- phylumData[,-1]})

# agregate counts at class level
combinedTaxa <- lapply(list_taxData,function(taxData) combineTaxa2(taxData[,-8],rank="class",confidence=0.8,returnFull=T))
classData <- lapply(seq_along(combinedTaxa),function(i) combCounts(combinedTaxa[[i]],counts(list_dds[[i]],normalize=T)))
classTaxa <- lapply(seq_along(combinedTaxa),function(i) combTaxa(combinedTaxa[[i]],list_taxData[[i]][,-8]))
classData <- lapply(seq_along(classData), function(i) aggregate(classData[[i]],by=list(taxaConfVec(classTaxa[[i]],conf=0.8,level=2)),sum))
classData <- lapply(classData,function(classData) {rownames(classData) <- classData[,1];classData <- classData[,-1]})


fl <- list(~site+loc_factor+condition*time*genotype_name,
					 ~loc_factor+condition*time*genotype_name,
					 ~loc_factor+condition*time*genotype_name,
					 ~loc_factor+time*genotype_name,
					 ~loc_factor+time*genotype_name,
					 ~loc_factor+time*genotype_name,
					 ~loc_factor+time*genotype_name)

myrda <- lapply(seq_along(list_dds),function(i) rda(t(phylumData[[i]]),formula=fl[[i]],data=colData(list_dds[[i]])))

species <- lapply(myrda,function(myrda) as.data.frame(scores(myrda,scaling="symmetric")$species))
species <- lapply(species,function(species) {species$phylum<-rownames(species);species})
# filter out phyla with less than 2% of total abundance
species <- lapply(seq_along(species), function(i )species[[i]][((rowSums(phylumData[[i]])/sum(phylumData[[i]]))*100)>=2,])
#centroids <- centroids[centroids$phylum %in% rhiz$phylum[rhiz$all>=4],]
sites   <- lapply(myrda,function(myrda) scores(myrda)$sites)

design <- c("site",rep("condition",2),rep("genotype_name",4))
titles <- c("Both orchards","Cider orchard","Dessert orchard","Cider tree station","Cider grass alley","Dessert tree station","Dessert grass alley")

exclude <- lapply(sites,function(sites) rownames(sites[!outliers::scores(sites, type="z", prob=0.9999)[,1]|outliers::scores(sites, type="z", prob=0.9999)[,2],]))

gl <- lapply(seq_along(sites),function(i)
						plotOrd(sites[[i]],list_dds[[i]]@colData,design=design[i],shape="time",axes=c(1,2),alpha=0.75,legend=T,continuous=F,exclude=exclude[i])+
						ggtitle(titles[i]) +
						geom_segment(data=species[[i]],aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
						geom_label(data=species[[i]],aes(label=phylum,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))

ggsave(paste(RHB,"RDA_phylum_biplot.pdf",sep="_"),plotOrd(sites,dds@colData,design="Cultivar",shape="Sus",axes=c(1,2),alpha=0.75,legend=T,continuous=F,ylims=c(-100,100),xlims=c(-100,100))  +
 # geom_point(data=species,aes(x=RDA1,y=RDA2,fill=phylum,shape=phylum), size=4,inherit.aes=F) +
	geom_segment(data=species,aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
	geom_label(data=species,aes(label=phylum,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))

ggsave(paste(RHB,"RDA_phylum_biplot_continuous.pdf",sep="_"),plotOrd(sites,dds@colData,design="Susceptibility",axes=c(1,2),alpha=0.75,legend=T,continuous=T,ylims=c(-100,100),xlims=c(-100,100))  +
 # geom_point(data=species,aes(x=RDA1,y=RDA2,fill=phylum,shape=phylum), size=4,inherit.aes=F) +
	geom_segment(data=species,aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
	geom_label(data=species,aes(label=phylum,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))


test <- rda(t(phylumData


pdf(paste(RHB,"PCA.pdf",sep="_"))
 plotOrd(d[[1]],colData(list_dds[[1]]),design="site",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Both orchards") +
 plotOrd(d[[2]],colData(list_dds[[2]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider orchard")
 plotOrd(d[[3]],colData(list_dds[[3]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Dessert orchard")
 plotOrd(d[[4]],colData(list_dds[[4]]),design="genotype_name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75) + ggtitle("Cider tree station")
 plotOrd(d[[5]],colData(list_dds[[5]]),design="genotype_name",shape="time"
