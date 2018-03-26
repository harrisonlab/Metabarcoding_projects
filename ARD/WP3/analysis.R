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
dds <- dds[,gsub("(^[A-Z][0-9]*)([A-Z])(.*)","\\2",rownames(colData(dds)))!="C"]

# There are only 3 (out of 900) missing samples - subsampling is a bit too extreme

# get number of samples per tree
# sample_numbers <- table(sub("[A-Z]$","",rownames(colData(dds))))

# collapse samples to mean
dds <- collapseReplicates2(dds,simple=T,groupby=sub("[A-Z]$","",rownames(colData(dds))))

# set the dds sizefactor to the number of samples 
# dds$sizeFactors <- as.vector(3/sample_numbers)

# recreate countData and colData
countData<- counts(dds,normalize=F)
colData <- as.data.frame(colData(dds))

# new dds object with the corrected data set
# dds <- DESeqDataSetFromMatrix(countData,colData,~1)

# calculate size factors - using geoMeans function (works better with this data set)
max(geoMeans(dds))/min(geoMeans(dds))
max(sizeFactors(estimateSizeFactors(dds)))/min(sizeFactors(estimateSizeFactors(dds)))
# sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
sizeFactors(dds) <-geoMeans(dds) 

#===============================================================================
#       subset data
#===============================================================================

# subset into orchards
dds_cider <- dds[,substr(rownames(colData(dds)),1,1)=="H"]
dds_dessert <- dds[,substr(rownames(colData(dds)),1,1)=="G"]

# subset cider orchard into trees and grass allyways
dds_cider_tree   <- dds_cider[,dds_cider$condition=="Y"]
dds_cider_grass  <- dds_cider[,dds_cider$condition=="N"]

# subset dessert orchard into trees and grass allyways
dds_dessert_tree   <- dds_dessert[,dds_dessert$condition=="Y"]
dds_dessert_grass  <- dds_dessert[,dds_dessert$condition=="N"]

# create a list of all dds objects to make it easier to apply methods
dds_list <- list(dds=dds,
		 dds_cider=dds_cider,
		 dds_dessert=dds_dessert,
		 dds_cider_tree=dds_cider_tree,
		 dds_cider_grass=dds_cider_grass,
		 dds_dessert_tree=dds_dessert_tree,
		 dds_dessert_grass=dds_dessert_grass)

#===============================================================================
#       Alpha diversity analysis - RUN BEFORE FILTERING OUT ANY LOW COUNT OTUS
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

# plot alpha diversity - all data
ggsave(paste(RHB,"Alpha_choa_all.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,4000,"Chao1")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",returnData=T)

# join diversity indices and metadata
dds$Samples <- rownames(colData(dds)) # or could use tibble/rownames construct in join by syntax)
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,as.data.frame(colData(dds))))

# remove specultation site (not necessary as I now have 
# all_alpha_ord <- all_alpha_ord[site!="Speculation",]

# perform anova for each index (this may need editing as the design will be unbalanced)
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
	setkey(all_alpha_ord,S.chao1)
	print("Chao1")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~condition*site,all_alpha_ord))
	setkey(all_alpha_ord,shannon)
	print("Shannon")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~condition*site,all_alpha_ord))
	setkey(all_alpha_ord,simpson)
	print("simpson")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~condition*site,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#===============================================================================

dds_list <- lapply(dds_list,function(o) {o[rowSums(counts(o, normalize=T))>4,]})

#===============================================================================
#       Beta diversity analysis
#===============================================================================

### PCA ###

# perform PC decomposition of DES objects
mypca_list <- lapply(dds_list,des_to_pca)		   

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d_list <- lapply(mypca_list,function(o) {t(data.frame(t(o$x)*o$percentVar))})

# plot the PCAs
pdf("PCA_plots.pdf",width=7,height=7)
# all data
plotOrd(d_list[[1]],colData(dds_list[[1]]),design="site",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
# orchard data
plotOrd(d_list[[2]],colData(dds_list[[2]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
plotOrd(d_list[[3]],colData(dds_list[[3]]),design="condition",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
# row data
plotOrd(d_list[[4]],colData(dds_list[[4]]),design="genotype.name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
plotOrd(d_list[[5]],colData(dds_list[[5]]),design="genotype.name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
plotOrd(d_list[[6]],colData(dds_list[[6]]),design="genotype.name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
plotOrd(d_list[[7]],colData(dds_list[[7]]),design="genotype.name",shape="time",pointSize=1.5,axes=c(1,2),alpha=0.75)
dev.off()
# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
	print("ANOVA")
	lapply(seq(1:4),function(x) {
		summary(aov(mypca$x[,x]~condition*site,colData(dds)))
	})
	print("PERMANOVA")
	lapply(seq(1:4),function(x) {
		summary(aovp(mypca$x[,x]~condition*site,colData(dds)))
	})
sink()

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

###	CCA ###

ord_cca <- ordinate(myphylo,method="CCA","samples",formula=~Treatment + Genotype + Treatment * Genotype + Condition(Block))

plot_ordination(myphylo, ord_cca, "samples", color="Treatment",shape="Genotype")

anova.cca(ord_cca)

### RDA ###

# transform data using vst
otu_table(myphylo) <-  otu_table(assay(varianceStabilizingTransformation(dds),taxa_are_rows=T)

# calculate rda1 (treatment + genotype)
ord_rda1 <- ordinate(myphylo,method="RDA","samples",formula=~Treatment + Genotype)

# calculate rda2 (treatment + genotype + interaction)
ord_rda2 <- ordinate(myphylo,method="RDA","samples",formula=~Treatment + Genotype + Treatment * Genotype)

# permutation anova of rda1 and rda 2
aov_rda1 <- anova.cca(ord_rda1,permuations=9999)
aov_rda2 <- anova.cca(ord_rda2,permuations=9999)

## partial RDA

# calculate rda3 removing block effect(treatment + genotype)
ord_rda3 <- ordinate(myphylo,method="RDA","samples",formula= ~Condition(Block) + Treatment + Genotype)

# calculate rda4 removing block effect(treatment + genotype + interaction)
ord_rda4 <- ordinate(myphylo,method="RDA","samples",formula= ~Condition(Block) + Treatment + Genotype + Treatment * Genotype)

# permutation anova of rda3 and rda 4
aov_rda3 <- anova.cca(ord_rda3,permuations=9999)
aov_rda4 <- anova.cca(ord_rda4,permuations=9999)

# plots

p1 <- plot_ordination(myphylo, ord_rda1, "samples", color="Treatment",shape="Genotype")
p2 <- plot_ordination(myphylo, ord_rda2, "samples", color="Treatment",shape="Genotype")
p3 <- plot_ordination(myphylo, ord_rda3, "samples", color="Treatment",shape="Genotype")
p4 <- plot_ordination(myphylo, ord_rda4, "samples", color="Treatment",shape="Genotype")

ggsave(paste(RHB,"RDA1.pdf",sep="_"),p1)
ggsave(paste(RHB,"RDA2.pdf",sep="_"),p2)
ggsave(paste(RHB,"RDA3.pdf",sep="_"),p3)
ggsave(paste(RHB,"RDA4.pdf",sep="_"),p4)

ggsave(paste(RHB,"RDA1_facet.pdf",sep="_"),p1+ facet_wrap(~Genotype, 2)+ geom_point(size=3.5,alpha=0.75))
ggsave(paste(RHB,"RDA2_facet.pdf",sep="_"),p2+ facet_wrap(~Genotype, 2)+ geom_point(size=3.5,alpha=0.75))
ggsave(paste(RHB,"RDA3_facet.pdf",sep="_"),p3+ facet_wrap(~Genotype, 2)+ geom_point(size=3.5,alpha=0.75))
ggsave(paste(RHB,"RDA4_facet.pdf",sep="_"),p4+ facet_wrap(~Genotype, 2)+ geom_point(size=3.5,alpha=0.75))

sink(paste(RHB,"RDA_permutation_anova",sep="_"))
  print(aov_rda1)
  print(aov_rda2)
  print(aov_rda3)
  print(aov_rda4)
sink()
