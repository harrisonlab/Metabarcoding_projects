#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
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
library(metafuncs)
library(BiocParallel)

register(MulticoreParam(12))
environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)

#===============================================================================
#       Load data (and clean it up)
#===============================================================================

ubiom_BAC <- loadData("BAC.otus_table.txt","colData","BAC.taxa","BAC.phy",RHB="BAC")
ubiom_BAC$countData <- ubiom_BAC$countData[,colnames(ubiom_BAC$countData)%in%rownames(ubiom_BAC$colData)]
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
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

# remove ungrafted samples (if not required in analysis)
dds <- dds[,dds$Genotype!="M9_ungrafted"]
dds$Genotype <- droplevels(dds$Genotype)

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
qp <- function(obj,name,colData,axes=c(1,2)) {
  ggsave(paste0(RHB,"_",name,".pdf"),plotOrd(obj,colData,design="Treatment",shape="Genotype",pointSize=1.5,alpha=0.75,axes=axes)+theme_classic_thin())
  ggsave(paste0(RHB,"_",name,"_facet.pdf"),plotOrd(obj,colData,design="Treatment",facet="Genotype",pointSize=1.5,alpha=0.75,axes=axes)+facet_wrap(~facet,3)+theme_facet_blank(angle=0))
  ggsave(paste0(RHB,"_",name,"_facet_bw.pdf"),plotOrd(obj,colData,shapes="Treatment",facet="Genotype",pointSize=1.5,alpha=0.75,axes=axes)+facet_wrap(~facet,3)+theme_facet_blank(angle=0))	
}
axes=c(1,2)
qp(d,"PCA",colData(dds),axes)


### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,as.data.frame(colData(dds))))

# add tree to phyloseq object
phy_tree(myphylo) <- nj(as.dist(phylipData))

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

theme_set(theme_bw())
p1 <- plot_ordination(myphylo, ordu, type="Samples", color="Treatment",shape="Genotype")
p1 + facet_wrap(~Genotype, 3)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData(dds),design="Block",xlabel="NMDS1",ylabel="NMDS2",pointSize=2),width=10,height=10)

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
 print("weighted")
 adonis(distance(myphylo,"unifrac",weighted=T)~Block + Treatment + Genotype + Treatment * Genotype,colData(dds),parallel=12,permutations=9999)
 print("unweighted")
 adonis(distance(myphylo,"unifrac",weighted=F)~Block + Treatment + Genotype + Treatment * Genotype,colData(dds),parallel=12,permutations=9999)
sink()

#===============================================================================
#      Population structure CCA/RDA
#===============================================================================

myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,as.data.frame(colData(dds))))
       
### CCA ###

ord_cca <- ordinate(myphylo,method="CCA","samples",formula=~Treatment + Genotype + Treatment * Genotype + Condition(Block))

plot_ordination(myphylo, ord_cca, "samples", color="Treatment",shape="Genotype")

anova.cca(ord_cca)

### RDA ###

# transform data using vst
otu_table(myphylo) <-  otu_table(assay(varianceStabilizingTransformation(dds)),taxa_are_rows=T)

# calculate rda1 (treatment + genotype)
ord_rda1 <- ordinate(myphylo,method="RDA","samples",formula=~Treatment + Genotype)

# calculate rda2 (treatment + genotype + interaction)
ord_rda2 <- ordinate(myphylo,method="RDA","samples",formula=~Treatment + Genotype + Treatment * Genotype)

# permutation anova of rda1 and rda 2
aov_rda1 <- anova.cca(ord_rda1,permuations=1000)
aov_rda2 <- anova.cca(ord_rda2,permuations=1000)

## partial RDA

# calculate rda3 removing block effect(treatment + genotype)
ord_rda3 <- ordinate(myphylo,method="RDA","samples",formula= ~Condition(Block) + Treatment + Genotype)

# calculate rda4 removing block effect(treatment + genotype + interaction)
ord_rda4 <- ordinate(myphylo,method="RDA","samples",formula= ~Condition(Block) + Treatment + Genotype + Treatment * Genotype)

# permutation anova of rda3 and rda 4
aov_rda3 <- anova.cca(ord_rda3,permuations=1000)
aov_rda4 <- anova.cca(ord_rda4,permuations=1000)

# sig
sink(paste(RHB,"RDA_permutation_anova",sep="_"))
 print(aov_rda1)
 print(aov_rda2)
 print(aov_rda3)
 print(aov_rda4)
sink()

### plots ###
	
#scores scaled by variation in each axes
sscores <- function(ord,axes=c(1,2)) {
	d <- scores(ord)$sites 
	eigvec = eigenvals(ord)
	fracvar = eigvec[axes]/sum(eigvec)
	percVar = fracvar * 100
	d <- t(t(d)*percVar)
	d
}

axes=c(1,2)
qp(sscores(ord_rda1),"RDA1",colData(dds),axes)
qp(sscores(ord_rda2),"RDA2",colData(dds),axes)
qp(sscores(ord_rda3),"RDA3",colData(dds),axes)
qp(sscores(ord_rda4),"RDA4",colData(dds),axes)


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
	
