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


ubiom_BAC <- loadData("BAC.otus_table.txt","colData","BAC.taxa","BAC.phy",RHB="BAC")
ubiom_FUN <- loadData("FUN.otus_table.txt","colData","FUN.taxa","FUN.phy",RHB="FUN")
ubiom_OO <- loadData("OO.otus_table.txt","colData2","OO.taxa","OO.phy",RHB="OO")
ubiom_NEM <- loadData("NEM.otus_table.txt","colData2","NEM.taxa","NEM.phy",RHB="NEM")

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
#       Create DEseq objects 
#===============================================================================

ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC,filter=expression(colSums(countData)>=1000))
ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN,filter=expression(colSums(countData)>=1000))
ubiom_OO$dds <- ubiom_to_des(ubiom_OO,filter=expression(colSums(countData)>=1000))
ubiom_NEM$dds <- ubiom_to_des(ubiom_NEM,filter=expression(colSums(countData)>=1000))

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# Recreate dds object and don't filter for low counts before running Alpha diversity

# alphaDiversity(ubiom_BAC,design="Treatment",model=~Genotype+Treatment+Genotype*Treatment,measures=c("Chao1", "Shannon", "Simpson","Observed")
# nice to do above, but time is limited so will stick to normal way of doing things (or list everything out)

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values 
ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData,design="Condition",colour="Location",measures=c("Chao1", "Shannon", "Simpson","Observed")))

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
phy_tree(myphylo) <- nj(as.dist(phylipData)) 

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData,shape="Condition",design="Location",continuous=T,xlabel="NMDS1",ylabel="NMDS2",alpha=0.75,pointSize=2))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
  adonis(distance(myphylo,"unifrac",weighted=T)~Pair+Condition,colData(dds),parallel=12,permutations=9999)
sink()

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

### CCA/RDA ###         
         
         
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
