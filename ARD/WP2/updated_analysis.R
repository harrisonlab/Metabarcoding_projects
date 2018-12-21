#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(data.table)
library(tidyverse)
library(Biostrings)
library(vegan)
library(lmPerm)
library(phyloseq)
library(ape)
#library(devtools)
#install_github("eastmallingresearch/Metabarcoding_pipeline/scripts") 
library(metafuncs)
library(BiocParallel)

register(MulticoreParam(12))
environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)

#===============================================================================
#       Load data (and clean it up)
#===============================================================================
# not interested in bacteria
# ubiom_BAC <- loadData("BAC.otus_table.txt","colData","BAC.taxa","BAC.phy",RHB="BAC")
# ubiom_BAC$countData <- ubiom_BAC$countData[,colnames(ubiom_BAC$countData)%in%rownames(ubiom_BAC$colData)]
ubiom_FUN <- loadData("FUN.otus_table.txt","colData","FUN.taxa","FUN.phy",RHB="FUN")
ubiom_FUN$countData <- ubiom_FUN$countData[,colnames(ubiom_FUN$countData)%in%rownames(ubiom_FUN$colData)]
ubiom_OO <- loadData("OO.otus_table.txt","colData","OO.taxa","OO.phy",RHB="OO")
rownames(ubiom_OO$colData) <- ubiom_OO$colData$Sample_ON
ubiom_NEM <- loadData("NEM.otus_table.txt","colData","NEM.taxa","NEM.phy",RHB="NEM")
rownames(ubiom_NEM$colData) <- ubiom_NEM$colData$Sample_ON
ubiom_NEM$countData <- ubiom_NEM$countData[,colnames(ubiom_NEM$countData)%in%rownames(ubiom_NEM$colData)]

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

# oomycetes
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("OO.taxa")
combinedTaxa <- combinedTaxa[c(1,3,5),]
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_OO$countData <- countData
ubiom_OO$taxData <- taxData

# Nematodes
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("NEM.taxa")
combinedTaxa <- combinedTaxa[1,]
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_NEM$countData <- countData
ubiom_NEM$taxData <- taxData

#===============================================================================
#       Create DEseq objects
#===============================================================================

ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN,filter=expression(colSums(countData)>=1000&colData$Block!="R"))
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC,filter=expression(colSums(countData)>=1000&colData$Block!="R"))
ubiom_OO$dds <- ubiom_to_des(ubiom_OO,filter=expression(colSums(countData)>=1000&colData$Block!="R"),calcFactors=geoMeans)
ubiom_NEM$dds <- ubiom_to_des(ubiom_NEM,filter=expression(colSums(countData)>=500&colData$Block!="R"))

#===============================================================================
#     Nematodes  Filter data
#===============================================================================

invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))

myfilter <- row.names(taxData[as.number(taxData$c_conf)>0.9 & as.number(taxData$o_conf)>0.9,])
ubiom_NEM$dds <- dds[rownames(dds)%in%myfilter,]

#===============================================================================
#     OOMYCETES  Filter data
#===============================================================================

invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))

### OOMYCETE FILTER to remove OTUs which are unlikely part of the correct kingdom (best to do this before Alpha diversity analysis)
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
ubiom_OO$dds <- dds[myfilter,]

#===============================================================================
#     Common pipeline
#===============================================================================

# attach objects (On of FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
#invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

# remove ungrafted samples (if not required in analysis)
dds <- dds[,dds$Genotype!="M9_ungrafted"]
dds$Genotype <- droplevels(dds$Genotype)

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
#       Phylum (or other taxa) plots
#===============================================================================

dd <- sumTaxa(list(data.frame(cbind(taxData[,1,drop=F],1)[,2,drop=F]),taxData,data.frame(all=1)),conf=0.9,proportional=T,taxon="phylum")
colnames(dd)[2] <- "all"
dd <- rbind(dd[dd$all>=1,],c("others",sum(dd[dd$all<1,2])))
md <- melt(dd,id=colnames(dd)[1])
#levels(md$variable)[1] <- "OTUs"
md$value <- as.numeric(md$value)
#md$Normalisation <- md$variable
md$phylum <- sub(".*_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)
g <- ggplot(md,aes(x=phylum,y=value))
g <- g + geom_bar(stat="identity",colour="black")
g <- g  + xlab("")+ ylab("")
scaleFUN<-function(x) sprintf("%.0f", x)
g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=NULL)
g <- g + guides(fill=guide_legend(ncol=1))
g <- g + theme_blank()
g <- g + theme(
	axis.text.x = element_text(angle = 45, hjust = 1,size=14),
	plot.margin=unit(c(0.2,0,0.2,1.5),"cm"),
	axis.line.y = element_line(colour = "black",size=1),
	axis.ticks.x=element_blank(),
	text=element_text(size=14),
	plot.title = element_text(hjust = -0.11),
	axis.title.y=element_text(size=(14-2)))
ggsave(paste0(RHB,"_OTU_frequency.pdf"),g,width=8,height=7)


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
#===============================================================================

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:5),function(x) summary(aov(mypca$x[,x]~Block + Treatment + Genotype + Treatment * Genotype,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:5),function(x) summary(aovp(mypca$x[,x]~Block + Treatment + Genotype + Treatment * Genotype,colData(dds))))
sink()


# plot  PCA
qp <- function(obj,name,colData) {
  ggsave(paste0(RHB,"_",name,".pdf"),plotOrd(obj,colData,design="Treatment",shape="Genotype",pointSize=1.5,alpha=0.75)+theme_classic_thin())
  ggsave(paste0(RHB,"_",name,"_facet.pdf"),plotOrd(obj,colData,design="Treatment",facet="Genotype",pointSize=1.5,alpha=0.75)+facet_wrap(~facet,3)+theme_facet_blank(angle=0))
  ggsave(paste0(RHB,"_",name,"_facet_bw.pdf"),plotOrd(obj,colData,facet="Treatment",shapes="Genotype",pointSize=1.5,alpha=0.75)+facet_wrap(~facet,3)+theme_facet_blank(angle=0))	
}
axes=c(1,2)
qp(d[,axes],"PCA",colData(dds))

#===============================================================================
#      Population structure CCA/RDA
#===============================================================================

myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,as.data.frame(colData(dds))))

### RDA ###

# transform data using vst
otu_table(myphylo) <-  otu_table(assay(varianceStabilizingTransformation(dds)),taxa_are_rows=T)
	
# no faf (consistent with others) method
sink(paste(RHB,"RDA_permutation_anova",sep="_"))
 ord_rda <- ordinate(myphylo,method="RDA","samples",formula= ~Block + Treatment + Genotype + Treatment * Genotype)		
 print(anova.cca(ord_rda,permuations=1000))
 print(anova.cca(ord_rda,permuations=1000,by="terms"))
sink()
sink(paste(RHB,"Partial_RDA_permutation_anova",sep="_"))
 ord_rda_partial <- ordinate(myphylo,method="RDA","samples",formula= ~Condition(Block) + Treatment + Genotype + Treatment * Genotype)
 print(anova.cca(ord_rda_partial,permuations=1000))
 print(anova.cca(ord_rda_partial,permuations=1000,by="terms"))
sink()

### plots ###
	
#scores scaled by variation in each axes
sscores <- function(ord,axes=c(1,2)) {
	d <- scores(ord,axes)$sites 
	eigvec = eigenvals(ord)
	fracvar = eigvec[axes]/sum(eigvec)
	percVar = fracvar * 100
	d <- t(t(d)*percVar)
	d
}

qp(sscores(ord_rda,c(2,3)),"RDA",colData(dds))
qp(sscores(ord_rda_partial),"Partial_RDA",colData(dds))

#===============================================================================
#       differential analysis
#===============================================================================

# p value for FDR cutoff
alpha <- 0.1

#### Treatment effect
design(dds) <- ~Block + Treatment
dds <- DESeq(dds,parallel=T)

res1 <- results(dds,alpha=alpha,parallel=T,contrast=c("Treatment","Nematicide","Control"))
res2 <- results(dds,alpha=alpha,parallel=T,contrast=c("Treatment","Nem_Fung","Control"))
res3 <- results(dds,alpha=alpha,parallel=T,contrast=c("Treatment","Nem_Oom","Control"))
res4 <- results(dds,alpha=alpha,parallel=T,contrast=c("Treatment","Nem_Oom_Fung","Control"))

summary(res1)
summary(res2)
summary(res3)
summary(res4)

#ggsave(paste(RHB,"MA_N.pdf",sep="_"),plot_ma(res1[,c(2,1,6)],legend=T))
#ggsave(paste(RHB,"MA_N_F.pdf",sep="_"),plot_ma(res2[,c(2,1,6)],legend=T))
#ggsave(paste(RHB,"MA_N_O.pdf",sep="_"),plot_ma(res3[,c(2,1,6)],legend=T))
#ggsave(paste(RHB,"MA_N_O_F.pdf",sep="_"),plot_ma(res4[,c(2,1,6)],legend=T))

test <- aggregate(t(counts(dds,normalize=T)),list(dds$Treatment),sum)	
# problem with incorrect reporting of FC and sig values when all data for an OTU in contrast is zero
# This issue is a known (mostly) none issue. Solution is to use LRT (See https://support.bioconductor.org/p/104803/)	
# sum(apply(test[,-1],2,function(x) {x[1]==0&product(x[-1])==0}))
# test[,c(T,apply(test[,-1],2,function(x) {x[1]==0&product(x[-1])==0}))]
# test[,c("Group.1","OTU248","OTU308","OTU367")]
# colnames(test[,test[1,]==0&test[2,]==0])

res1$trustworthy=1;res2$trustworthy=1;res3$trustworthy=1;res4$trustworthy=1;
res1[colnames(test[,test[1,]==0&test[2,]==0]),]$trustworthy=0
res2[colnames(test[,test[1,]==0&test[3,]==0]),]$trustworthy=0
res3[colnames(test[,test[1,]==0&test[4,]==0]),]$trustworthy=0
res4[colnames(test[,test[1,]==0&test[5,]==0]),]$trustworthy=0	

res1 <- as.data.table(as.data.frame(res1),keep.rownames="OTU")
res2 <- as.data.table(as.data.frame(res2),keep.rownames="OTU")
res3 <- as.data.table(as.data.frame(res3),keep.rownames="OTU")
res4 <- as.data.table(as.data.frame(res4),keep.rownames="OTU")	
	
taxDT <- as.data.table(taxData,keep.rownames="OTU")	
	
## LRT TESTS ##
# dds <- estimateDispersions(dds) # if not already calculated 
	
full <- model.matrix(design(dds), colData(dds))
reduced <- subset(full,select=-TreatmentNematicide)
dds <- nbinomLRT(dds, full=full, reduced=reduced)
res <- as.data.table(as.data.frame(results(dds)),keep.rownames="OTU")
fwrite(res[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_LRT.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
res1[,LRT_padj:=res1[res[,c(1,7)],on="OTU",nomatch=0][,9]]
	
reduced <- subset(full,select=-TreatmentNem_Fung)
dds <- nbinomLRT(dds, full=full, reduced=reduced)
res <- as.data.table(as.data.frame(results(dds)),keep.rownames="OTU")
fwrite(res[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_Fung_LRT.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
res2[,LRT_padj:=res2[res[,c(1,7)],on="OTU",nomatch=0][,9]]

reduced <- subset(full,select=-TreatmentNem_Oom)
dds <- nbinomLRT(dds, full=full, reduced=reduced)
res <- as.data.table(as.data.frame(results(dds)),keep.rownames="OTU")
fwrite(res[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_Oom_LRT.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
res3[,LRT_padj:=res3[res[,c(1,7)],on="OTU",nomatch=0][,9]]

reduced <- subset(full,select=-TreatmentNem_Oom_Fung)
dds <- nbinomLRT(dds, full=full, reduced=reduced)
res <- as.data.table(as.data.frame(results(dds)),keep.rownames="OTU")
fwrite(res[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_Oom_Fung_LRT.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
res4[,LRT_padj:=res4[res[,c(1,7)],on="OTU",nomatch=0][,9]]

fwrite(res1[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nematicide_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
fwrite(res2[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_Fung_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
fwrite(res3[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_Oom_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
fwrite(res4[taxDT,nomatch=0,on="OTU"],paste(RHB,"Nem_Oom_Fung_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

### treatment effect per genotype ###

# the full model (looking for effect of treatment per genotype) - o.k. this is easier if I combine Treatment and Genotype
#full_design <- ~Block + Genotype + Treatment + Genotype * Treatment

# add full model to dds object
#design(dds) <- full_design

# add grouping vector
dds$Group <- as.factor(paste(dds$Genotype,dds$Treatment,sep="."))

# add design to dds object
design(dds) <- ~Block + Group

# calculate fit using precalculated dispersions
#dds <- nbinomWaldTest(dds)
dds <- DESeq(dds,parallel=T)

# Treatment effects 
# contrast <- c("Group","G41.Nem_Fung","G41.Control")
# contrast <- c("Group","M9.Nem_Fung","M9.Control")

# calculate results (will output 4 x 5 matrix: cols genotype, rows treatment)
res <- 
sapply(levels(dds$Genotype),function(x) {
 sapply(levels(dds$Treatment)[-1],function(y) {
  treatment <- paste(x,y,sep=".")
  control <- paste(x,"Control",sep=".")
  results(dds,alpha=alpha,parallel=T,contrast=c("Group",treatment,control))  
 })  
})

# output results to files
sapply(seq_along(res),function(i) write.table(data.table(inner_join(data.table(OTU=rownames(res[[i]]),as.data.frame(res[[i]])),data.table(OTU=rownames(taxData),taxData))),paste(RHB, sub(".*Group ","",res[[i]]@elementMetadata$description[[2]]),"txt",sep="."),quote=F,sep="\t",na="",row.names=F))


# output sig fasta
# writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa"))
	
