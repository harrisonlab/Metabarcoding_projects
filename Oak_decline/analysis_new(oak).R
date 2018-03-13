#===============================================================================
#       Load libraries
#===============================================================================
library(DESeq2)
library(phyloseq)
library(BiocParallel)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
register(MulticoreParam(12))
library(lmPerm)
library(ape)
library(vegan)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
environment(plot_alpha) <-environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)

# library(future)
# plan(multiprocess)

#===============================================================================
#       Load data 
#===============================================================================

ubiom_BAC <- loadData("16S.otu_table.txt","colData","16S.taxa","16S.phy",RHB="BAC")
ubiom_FUN <- loadData("ITS.otu_table.txt","colData","ITS.taxa","ITS.phy",RHB="FUN")

#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("ITS.taxa")
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData

#===============================================================================
#       Create DEseq objects (and remove none oak samples)
#===============================================================================

ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN,filter=expression(colData$OAK=="Quercus"),calcFactors=geoMeans)
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC,filter=expression(colData$OAK=="Quercus"),calcFactors=geoMeans)

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects (FUN, BAC,OO or NEM)
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Subsample data
#===============================================================================

# get number of samples per tree
sample_numbers <- table(paste0(dds$condition,dds$site,dds$tree))

# collapse (sum) samples
dds <- collapseReplicates(dds,groupby=paste0(dds$condition,dds$site,dds$tree))

# set sequence depth to sample
depth <- colSums(counts(dds))*(1/sample_numbers)

# set random seed  value to get repeatable results
fixed_seed <- sum(utf8ToInt("Xiangming Xu")) 

# sub sample counts 
sub_counts1 <- subsample(counts(dds),depth,replace=T,fixed_seed)
sub_counts2 <- subsample(counts(dds),depth,replace=F,fixed_seed) # sampling without replacement is probably a better scheme (though there's very little difference in the final results)

# create new deseq object
dds <- DESeqDataSetFromMatrix(sub_counts2,colData(dds),~1)

# calculate size factors - using geoMeans function (works better with this data set)
max(geoMeans(dds))/min(geoMeans(dds))
max(sizeFactors(estimateSizeFactors(dds)))/min(sizeFactors(estimateSizeFactors(dds)))
# sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
sizeFactors(dds) <-geoMeans(dds) 
# calcNormFactors(counts(dds),method="RLE",lib.size=(prop.table(colSums(counts(dds)))))

#===============================================================================
#       Alpha diversity analysis - RUN BEFORE FILTERING OUT ANY LOW COUNT OTUS
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values
ggsave(paste(RHB,"Alpha_all.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,600,"Chao1")))

### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",returnData=T)

# join diversity indices and metadata
dds$Samples <- rownames(colData(dds)) # or could use tibble/rownames construct in join by syntax)
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,as.data.frame(colData(dds))))

# remove specultation site..
all_alpha_ord <- all_alpha_ord[site!="Speculation",]

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

# chestnuts wood paired samples
all_alpha_ord <- all_alpha_ord[site=="Chestnuts",]
all_alpha_ord <- all_alpha_ord[c(tree[all_alpha_ord$condition=="Healthy" & all_alpha_ord$site=="Chestnuts"] %in% all_alpha_ord$tree[all_alpha_ord$condition=="Symptom" & all_alpha_ord$site=="Chestnuts"],all_alpha_ord$tree[all_alpha_ord$condition=="Symptom" & all_alpha_ord$site=="Chestnuts"] %in% all_alpha_ord$tree[all_alpha_ord$condition=="Healthy" & all_alpha_ord$site=="Chestnuts"])]
sink(paste(RHB,"chestnut_ALPHA_stats.txt",sep="_"))
	setkey(all_alpha_ord,S.chao1)
	print("Chao1")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~tree+condition,all_alpha_ord))
	setkey(all_alpha_ord,shannon)
	print("Shannon")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~tree+condition,all_alpha_ord))
	setkey(all_alpha_ord,simpson)
	print("simpson")
	summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~tree+condition,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#===============================================================================

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity analysis
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds[,dds$site!="Speculation"])

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
g <- plotOrd(d,colData(dds[,dds$site!="Speculation"]),design="site",shape="condition",pointSize=1.5,axes=c(1,2),alpha=0.75)
ggsave(paste(RHB,"PCA.pdf",sep="_"),g)

# ANOVA
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
	print("ANOVA")
	lapply(seq(1:4),function(x) {
		summary(aov(mypca$x[,x]~condition*site,colData(dds[,dds$site!="Speculation"])))
	})
	print("PERMANOVA")
	lapply(seq(1:4),function(x) {
		summary(aovp(mypca$x[,x]~condition*site,colData(dds[,dds$site!="Speculation"])))
	})
sink()

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds[,dds$site!="Speculation"],normalize=T),taxData,as.data.frame(colData(dds[,dds$site!="Speculation"]))))

# add tree to phyloseq object
phy_tree(myphylo) <- nj(as.dist(phylipData))

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData(dds[,dds$site!="Speculation"]),design="site",shape="condition",xlabel="NMDS1",ylabel="NMDS2",pointSize=1.5,axes=c(1,2),alpha=0.75))

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
	print("weighted")
	adonis(distance(myphylo,"unifrac",weighted=T)~condition*site,colData(dds[,dds$site!="Speculation"]),parallel=12,permutations=9999)
	print("unweighted")
	adonis(distance(myphylo,"unifrac",weighted=F)~condition*site,colData(dds[,dds$site!="Speculation"]),parallel=12,permutations=9999)

sink()

## including speculation site ##
mypca <- des_to_pca(dds)
d <-t(data.frame(t(mypca$x)*mypca$percentVar))
ggsave(paste(RHB,"ALL_SITES_PCA.pdf",sep="_"),plotOrd(d,colData(dds),design="site",shape="condition",pointSize=1.5,axes=c(1,2),alpha=0.75))
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,as.data.frame(colData(dds))))
phy_tree(myphylo) <- nj(as.dist(phylipData))
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
ggsave(paste(RHB,"Unifrac_ALL_SITES_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData(dds),design="site",shape="condition",xlabel="NMDS1",ylabel="NMDS2",pointSize=1.5,axes=c(1,2),alpha=0.75))

## just chestnuts ##
dds_nuts <- dds[,dds$site=="Chestnuts"]
dds_nuts <- dds_nuts[,c(dds$tree[dds$condition=="Healthy" & dds$site=="Chestnuts"] %in% dds$tree[dds$condition=="Symptom" & dds$site=="Chestnuts"],dds$tree[dds$condition=="Symptom" & dds$site=="Chestnuts"] %in% dds$tree[dds$condition=="Healthy" & dds$site=="Chestnuts"])]
mypca <- des_to_pca(dds_nuts)
sink(paste(RHB,"PCA_chestnuts_ANOVA.txt",sep="_"))
	print("ANOVA")
	lapply(seq(1:4),function(x) {
		summary(aov(mypca$x[,x]~tree+condition,colData(dds_nuts)))
	})
	print("PERMANOVA")
	lapply(seq(1:4),function(x) {
		summary(aovp(mypca$x[,x]~tree+condition,colData(dds_nuts)))
	})
sink()
myphylo <- ubiom_to_phylo(list(counts(dds_nuts,normalize=T),taxData,as.data.frame(colData(dds_nuts))))
phy_tree(myphylo) <- nj(as.dist(phylipData))
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
sink(paste(RHB,"chestnuts_PERMANOVA_unifrac.txt",sep="_"))
	print("weighted")
	adonis(distance(myphylo,"unifrac",weighted=T)~tree+condition,colData(dds_nuts),parallel=12,permutations=9999)
	print("unweighted")
	adonis(distance(myphylo,"unifrac",weighted=F)~tree+condition,colData(dds_nuts),parallel=12,permutations=9999)

sink()

#===============================================================================
#       differential analysis
#===============================================================================

# p value for FDR cutoff
alpha <- 0.5

# chestnuts wood paired samples
dds_nuts <- dds[,dds$site=="Chestnuts"]
dds_nuts <- dds_nuts[,c(dds$tree[dds$condition=="Healthy" & dds$site=="Chestnuts"] %in% dds$tree[dds$condition=="Symptom" & dds$site=="Chestnuts"],dds$tree[dds$condition=="Symptom" & dds$site=="Chestnuts"] %in% dds$tree[dds$condition=="Healthy" & dds$site=="Chestnuts"])]

# drop unused levels from colData
colData(dds_nuts) <- droplevels(colData(dds_nuts))

# the model
full_design <- ~tree+condition

# add model to dds object
design(dds_nuts) <- full_design

# calculate fit
dds_nuts <- DESeq(dds_nuts,parallel=T)	       

# calculate results for default contrast (S vs H)
res <- results(dds_nuts,alpha=alpha,parallel=T)

# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
write.table(res.merge, paste(RHB,"chestnut_paired_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)


# chestnut + bigwood samples
dds2 <- dds[,dds$site!="Speculation"]
# drop unused levels from colData
colData(dds2) <- droplevels(colData(dds2))
# design
full_design <- ~site+condition
# add model
design(dds2) <- full_design
# calculate fit
dds2 <- DESeq(dds2,parallel=T)	       
# calculate results
res <- results(dds2,alpha=alpha,parallel=T)
# merge results with taxonomy data
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
#output results
write.table(res.merge, paste(RHB,"condition_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)


# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa"))

## MA plots
plot_ma(res[,c(2,1,6)])

######### END OF UPDATED STUFF ########


