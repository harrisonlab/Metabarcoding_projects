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
#assignInNamespace("plot_ordination",value=plot_ordination,ns="phyloseq")

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

# get unifrac dist
phylipData <- fread.phylip("BAC.phy")

njtree <- nj(as.dist(phylipData))

# save data into a list
ubiom_BAC <- list(
  countData=countData,
  colData=colData,
  taxData=taxData,
  phylipData=phylipData,
  njtree=njtree,
  RHB="BAC"
)

# Fungi all in one call
ubiom_FUN <- list(
  countData=read.table("FUN.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  phylipData=fread.phylip("FUN.phy"),
  RHB="FUN"
)
ubiom_FUN$njtree <- nj(as.dist(ubiom_FUN$phylipData))

# Oomycetes
ubiom_OO <- list(
  countData=read.table("OO.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData2",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("OO.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  phylipData=fread.phylip("OO.phy"),
  RHB="OO"
)
ubiom_OO$njtree <- nj(as.dist(ubiom_OO$phylipData))
rownames(ubiom_OO$colData) <- paste0("X",gsub("_","\\.",ubiom_OO$colData$name),"_",sub("D.*","",rownames(ubiom_OO$colData)))
rownames(ubiom_OO$colData) <- sub("XG","G",rownames(ubiom_OO$colData))

# Nematodes
ubiom_NEM <- list(
  countData=read.table("NEM.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData2",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("NEM.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  phylipData=fread.phylip("NEM.phy"),
  RHB="NEM"
)
ubiom_NEM$njtree <- nj(as.dist(ubiom_NEM$phylipData))
rownames(ubiom_NEM$colData) <- paste0("X",gsub("_","\\.",ubiom_NEM$colData$name),"_",sub("D.*","",rownames(ubiom_NEM$colData)))
rownames(ubiom_NEM$colData) <- sub("XG","G",rownames(ubiom_NEM$colData))

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

# Nematodes
# oomycetes

#===============================================================================
#      ****FUNGI****
#===============================================================================

invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# remove low count and control samples
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"

# remove Pair of any sample with a low count
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)

# apply filter
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

#===============================================================================
#      **Normalised**
#===============================================================================

sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,1000,"Chao1")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"PCA_loc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))
                            
ggsave(paste(RHB,"PCA_Original.pdf",sep="_"),plotOrd(d,colData,design="Condition",continuous=F,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))
                          
g <- plotOrd(d,colData,design="Condition",continuous=F,axes=c(1,3),plot="Label",labelSize=2.5,cbPalette=T,label="Pair",legend="bottom")
g$layers[[1]] <- NULL
g  <- g + geom_point(size = 0, stroke = 0)  # OR  geom_point(shape = "") +
g  <- g + geom_label(show.legend = FALSE,size=2.5)
g  <- g + guides(colour = guide_legend(override.aes = list(size = 5, shape = c(utf8ToInt("H"), utf8ToInt("S")))))
g  <- g + scale_colour_manual(name = "Condition", breaks = c("H","S"), labels = c("",""),values=c("#000000", "#E69F00"))
ggsave(paste(RHB,"PCA_NEW_1vs3.pdf",sep="_"),g)

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_without_paired_var.pdf",sep="_"),plotOrd(d,colData,design="Condition",continuous=F,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
 adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"BRAY_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top7phyla = names(sort(phylum.sum, TRUE))[1:7]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top7phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

# calculate ordination
Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

# plot with plot_ordination
theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"NMDS_taxa_by_Condition.pdf",sep="_"))
 plot_ordination(H, Hord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
 plot_ordination(S, Sord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
dev.off()

# theme_set(theme_classic_thin())
# plot_ordination(myphylo, ordu, type="samples", shape="Condition", color="Location",continuous=T)

# PCoA
# ordu = ordinate(myphylo, "PCoA", "unifrac", weighted=TRUE)
# d <-t(data.frame(t(ordu$vectors)*ordu$values$Relative_eig))
# plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PCoA1",ylabel="PCoA2")

#===============================================================================
#       differential analysis
#===============================================================================

# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power
# but, no point keeping OTUs with 0 count
dds<-dds[ rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

# the full model
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[ res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa"))

#==============================================================================
#       **qPCR**
#===============================================================================

dds<-DESeqDataSetFromMatrix(countData,colData,design)

# Correction from aboslute quantification
sizeFactors(dds) <- 1/colData$funq

# Correction from aboslute quantification v2
# sizeFactors(dds) <- sizeFactors(dds)/colData$funq

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"qPCR_Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
 sink(paste(RHB,"qPCR_ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"qPCR_OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"qPCR_PCA.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"qPCR_PCA_loc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2,ylims=c(-2,4)))

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"qPCR_PCA_deloc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"qPCR_PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"qPCR_Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"qPCR_PERMANOVA_unifrac.txt",sep="_"))
 adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"qPC_BRAY_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top7phyla = names(sort(phylum.sum, TRUE))[1:7]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top7phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"qPCR_NMDS_taxa_by_Condition.pdf",sep="_"))
 plot_ordination(H, Hord, type="taxa", color="phylum",ylims=c(-0.8,1.2),xlims=c(-0.8,1.2))+ facet_wrap(~phylum, 3)
 plot_ordination(S, Sord, type="taxa", color="phylum",ylims=c(-0.8,1.2),xlims=c(-0.8,1.2))+ facet_wrap(~phylum, 3)
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================

# filter for low counts - this can affect the FD probability and DESeq2 does apply its own filtering for genes/otus with no power
# but, no point keeping OTUs with 0 count
dds<-dds[ rowSums(counts(dds,normalize=T))>0,]

# p value for FDR cutoff
alpha <- 0.1

# the full model
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"qPCR_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[ res.merge[ padj<=0.05]$OTU],paste0(RHB,".qPCR_sig.fa"))

#===============================================================================
#      ****BACTERIA****
#===============================================================================

invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# remove low count and control samples
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"

# remove Pair of any sample with a low count
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)

# apply filter
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

#===============================================================================
#      ****Normalised****
#===============================================================================

sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,2000,"Chao1")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
 sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

### read accumulation filter
# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"PCA_loc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_deloc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
 adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"BRAY_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top12phyla = names(sort(phylum.sum, TRUE))[1:12]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top12phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

# calculate ordination
Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

# plot with plot_ordination
theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"NMDS_taxa_by_Condition.pdf",sep="_"))
plot_ordination(H, Hord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
plot_ordination(S, Sord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
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
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa"))

#==============================================================================
#       ****qPCR****
#===============================================================================

dds<-DESeqDataSetFromMatrix(countData,colData,design)

# Correction from aboslute quantification
sizeFactors(dds) <- 1/colData$bacq

# Correction from aboslute quantification v2
# sizeFactors(dds) <- sizeFactors(dds)/colData$bacq

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"qPCR_Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,2000,"Chao1")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
 sink(paste(RHB,"qPCR_ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"qPCR_OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"qPCR_PCA.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"qPCR_PCA_loc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"qPCR_PCA_deloc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"qPCR_PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"qPCR_Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"qPCR_BRAY_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"qPCR_PERMANOVA_unifrac.txt",sep="_"))
  adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top12phyla = names(sort(phylum.sum, TRUE))[1:12]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top12phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"qPCR_NMDS_taxa_by_Condition.pdf",sep="_"))
  plot_ordination(H, Hord, type="taxa", color="phylum",ylims=c(-0.8,1.2),xlims=c(-0.8,1.2))+ facet_wrap(~phylum, 3)
  plot_ordination(S, Sord, type="taxa", color="phylum",ylims=c(-0.8,1.2),xlims=c(-0.8,1.2))+ facet_wrap(~phylum, 3)
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
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"qPCR_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".qPCR_sig.fa"))

#===============================================================================
#      ****OOMYCETES****
#===============================================================================

invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# remove low count and control samples
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"

# remove Pair of any sample with a low count
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)

# apply filter
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

### filter to remove OTUs which are unlikely part of the correct kingdom (best to do this before Alpha diversity analysis)
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"PCA_loc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_deloc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
 adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"BRAY_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top7phyla = names(sort(phylum.sum, TRUE))[1:7]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top7phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

# calculate ordination
Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

# plot with plot_ordination
theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"NMDS_taxa_by_Condition.pdf",sep="_"))
 plot_ordination(H, Hord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
 plot_ordination(S, Sord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
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
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa"))

#==============================================================================
#       **qPCR**
#===============================================================================

dds<-DESeqDataSetFromMatrix(countData,colData,design)

# Correction from aboslute quantification of fungal ITS 
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))/left_join(colData,ubiom_FUN$colData)$funq

### filter to remove OTUs which are unlikely part of the correct kingdom (best to do this before Alpha diversity analysis)
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"Alpha_qPCR.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
sink(paste(RHB,"ALPHA_stats_qPCR.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts_qPCR.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"PCA_qPCR.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"PCA_loc_qPCR.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_deloc_qPCR.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"PCA_ANOVA_qPCR.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS_qPCR.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac_qPCR.txt",sep="_"))
 adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"BRAY_NMDS_qPCR.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top7phyla = names(sort(phylum.sum, TRUE))[1:7]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top7phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

# calculate ordination
Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

# plot with plot_ordination
theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"NMDS_taxa_by_Condition_qPCR.pdf",sep="_"))
 plot_ordination(H, Hord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
 plot_ordination(S, Sord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
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
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff_qPCR.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".sig_qPCR.fa"))


#===============================================================================
#      ****NEMATODE****
#===============================================================================

invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# remove low count and control samples
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"

# remove Pair of any sample with a low count
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)

# apply filter
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]

# simple Deseq design
design<-~1

#create DES object
# colnames(countData) <- row.names(colData)
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

### filter to remove OTUs which are unlikely part of the correct kingdom (best to do this before Alpha diversity analysis)
myfilter <- row.names(taxData[as.number(taxData$c_conf)>0.9 & as.number(taxData$o_conf)>0.9,])
dds <- dds[rownames(dds)%in%myfilter,]

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Condition",cbPalette=T,legend="hidden",measures=c("Chao1", "Shannon", "Simpson","Observed")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData,design="Condition",returnData=T)

# add column names as row to metadata (or use tribble)
colData$samples <- rownames(colData)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(inner_join(all_alpha_ord,colData,by=c("Samples"="samples")))

# perform anova for each index
colData$Pair<-as.factor(colData$Pair)
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
 setkey(all_alpha_ord,S.chao1)
 print("Chao1")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,shannon)
 print("Shannon")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,simpson)
 print("simpson")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Pair+Condition,all_alpha_ord))
 setkey(all_alpha_ord,S.ACE)
 print("ACE")
 summary(aovp(as.numeric(as.factor(all_alpha_ord$S.ACE))~Pair+Condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

### read accumulation filter
# plot cummulative reads (will also produce a data table "dtt" in the global environment)
ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Add spatial information as a numeric and plot
colData$Location<-as.number(colData$Pair)

# plot the PCA
pdf(paste(RHB,"PCA.pdf",sep="_"))
 plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
 plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

ggsave(paste(RHB,"PCA_loc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2",alpha=0.75,pointSize=2))

### remove spatial information (this uses the factor "Pair" not the numeric "Location") and plot
pc.res <- resid(aov(mypca$x~colData$Pair,colData))
d <- t(data.frame(t(pc.res*mypca$percentVar)))
ggsave(paste(RHB,"PCA_deloc.pdf",sep="_"),plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2"))

# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
 print("ANOVA")
 lapply(seq(1:3),function(x) summary(aov(mypca$x[,x]~Pair+Condition,colData(dds))))
 print("PERMANOVA")
 lapply(seq(1:3),function(x) summary(aovp(mypca$x[,x]~Pair+Condition,colData(dds))))
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
 adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

# calculate NMDS ordination with bray-curtis distance matrix     
ordu = ordinate(myphylo, "NMDS", "bray",stratmax=50) 

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"BRAY_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# phyla plots - keep top 7 phyla only
phylum.sum = tapply(taxa_sums(myphylo), tax_table(myphylo)[, "phylum"], sum, na.rm=TRUE)
top7phyla = names(sort(phylum.sum, TRUE))[1:7]
myphylo_slim = prune_taxa((tax_table(myphylo)[, "phylum"] %in% top7phyla), myphylo)

# split phylo into H and S
H <- prune_samples(colData$Condition=="H",myphylo_slim)
S <- prune_samples(colData$Condition=="S",myphylo_slim)

# calculate ordination
Hord <- ordinate(H, "NMDS", "bray")
Sord <- ordinate(S, "NMDS", "bray")

# plot with plot_ordination
theme_set(theme_facet_blank(angle=0,vjust=0,hjust=0.5))
pdf(paste(RHB,"NMDS_taxa_by_Condition.pdf",sep="_"))
 plot_ordination(H, Hord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
 plot_ordination(S, Sord, type="taxa", color="phylum")+ facet_wrap(~phylum, 3)
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
full_design <- ~Pair + Condition

# add full model to dds object
design(dds) <- full_design

# calculate fit
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H)
res <- results(dds,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa"))
