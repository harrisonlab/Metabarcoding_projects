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

#===============================================================================
#       network analysis
#===============================================================================

library(cooccur)
library(future) 
plan(multiprocess) 


myfiltbiom <- prune_samples(sample_data(biom16)$site=="Chestnuts",biom16)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 


cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1 


fht16 <- future({cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)})
fst16 <- future({cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)})
#fs16 <- future({cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T))}

CHcoHmodel16 <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)
CHcoSmodel16 <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)
 
 
CHcoHmodel16$results$padj <- p.adjust(apply(CHcoHmodel16$results[,8:9],1, min),"BH") 
CHcoSmodel16$results$padj <- p.adjust(apply(CHcoSmodel16$results[,8:9],1, min),"BH") 
#CHcomodel$16results$padj <- p.adjust(apply(CHcomodel$results[,8:9],1, min),"BH") 


myfiltbiom <- prune_samples(sample_data(biom16)$site=="Bigwood",biom16) 
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom))) 
 
cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Control",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 

cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1  

BWfht16 <- future({cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)})
BWfst16 <- future({cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)})

BWcoHmodel16 <- value(BWfht16) #  cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T) #value(BWfth1) 
BWcoSmodel16 <- value(BWfst16) #  cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T) #value(BWfts1) 
#BWcomodel16 <-  cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T) 

BWcoHmodel16$results$padj <- p.adjust(apply(BWcoHmodel16$results[,8:9],1, min),"BH") 
BWcoSmodel16$results$padj <- p.adjust(apply(BWcoSmodel16$results[,8:9],1, min),"BH") 

##Speculation

myfiltbiom <- prune_samples(sample_data(biomITS)$site=="Speculation",biomITS)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))
cotable <- cotable[rowSums(cotable)>5,colSums(cotable)>5]

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 

cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1  


SPcoHmodel <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)
SPcoSmodel <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)


SPcoHmodel $results$padj <- p.adjust(apply(SPcoHmodel $results[,8:9],1, min),"BH") 
SPcoSmodel $results$padj <- p.adjust(apply(SPcoSmodel $results[,8:9],1, min),"BH") 


cotable[cotable>0] <- 1 

SPcomodel <- cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T)
SPcomodel$results$padj <- p.adjust(apply(SPcomodel$results[,8:9],1, min),"BH") 


myfiltbiom <- prune_samples(sample_data(biom16)$site=="Speculation",biom16)
cotable <- as.data.frame(as.matrix(otu_table(myfiltbiom)))

cotable_h <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Healthy",myfiltbiom)))] 
cotable_s <- cotable[,row.names(sample_data(prune_samples(sample_data(myfiltbiom)$condition=="Symptom",myfiltbiom)))] 

cotable_s <- cotable_s[rowSums(cotable_s)>5,colSums(cotable_s)>5]
cotable_h <- cotable_h[rowSums(cotable_h)>5,colSums(cotable_h)>5]
cotable_s[cotable_s>0] <- 1 
cotable_h[cotable_h>0] <- 1  


SPcoHmodel16 <- cooccur2(cotable_h,type = "spp_site",spp_names = T,thresh = T)
SPcoSmodel16 <- cooccur2(cotable_s,type = "spp_site",spp_names = T,thresh = T)

SPcoHmodel16 $results$padj <- p.adjust(apply(SPcoHmodel16 $results[,8:9],1, min),"BH") 
SPcoSmodel16 $results$padj <- p.adjust(apply(SPcoSmodel16 $results[,8:9],1, min),"BH") 

cotable <- cotable[rowSums(cotable)>5,colSums(cotable)>5]
cotable[cotable>0] <- 1 

SPcomodel16 <- cooccur2(cotable,type = "spp_site",spp_names = T,thresh = T)
SPcomodel16$results$padj <- p.adjust(apply(SPcomodel16$results[,8:9],1, min),"BH") 

## ITS
summary(CHcoHmodel)
summary(CHcoSmodel) 
summary(BWcoHmodel)
summary(BWcoSmodel)
summary(SPcomodel)

nrow(CHcoHmodel$results[CHcoHmodel$results$padj<=0.1&CHcoHmodel$results$p_gt<=0.1,])
nrow(CHcoHmodel$results[CHcoHmodel$results$padj<=0.1&CHcoHmodel$results$p_lt<=0.1,])
nrow(CHcoSmodel$results[CHcoSmodel$results$padj<=0.1&CHcoSmodel$results$p_gt<=0.1,])
nrow(CHcoSmodel$results[CHcoSmodel$results$padj<=0.1&CHcoSmodel$results$p_lt<=0.1,])
nrow(BWcoHmodel$results[BWcoHmodel$results$padj<=0.1&BWcoHmodel$results$p_gt<=0.1,])
nrow(BWcoHmodel$results[BWcoHmodel$results$padj<=0.1&BWcoHmodel$results$p_lt<=0.1,])
nrow(BWcoSmodel$results[BWcoSmodel$results$padj<=0.1&BWcoSmodel$results$p_gt<=0.1,])
nrow(BWcoSmodel$results[BWcoSmodel$results$padj<=0.1&BWcoSmodel$results$p_lt<=0.1,])
nrow(SPcoHmodel$results[SPcoHmodel$results$padj<=0.1&SPcoHmodel$results$p_gt<=0.1,])
nrow(SPcoHmodel$results[SPcoHmodel$results$padj<=0.1&SPcoHmodel$results$p_lt<=0.1,])
nrow(SPcoSmodel$results[SPcoSmodel$results$padj<=0.1&SPcoSmodel$results$p_gt<=0.1,])
nrow(SPcoSmodel$results[SPcoSmodel$results$padj<=0.1&SPcoSmodel$results$p_lt<=0.1,])


pdf("obs.exp.fungi.pdf")
obs.v.exp(CHcoHmodel)
obs.v.exp(CHcoHmodel)
obs.v.exp(BWcoHmodel)
obs.v.exp(BWcoSmodel)
obs.v.exp(SPcoHmodel)
obs.v.exp(SPcoSmodel)
dev.off()


## 16S

summary(CHcoHmodel16)
summary(CHcoSmodel16)
summary(BWcoHmodel16)
summary(BWcoSmodel16)
summary(SPcomodel16)

nrow(CHcoHmodel16$results[CHcoHmodel16$results$padj<=0.1&CHcoHmodel16$results$p_gt<=0.1,])
nrow(CHcoHmodel16$results[CHcoHmodel16$results$padj<=0.1&CHcoHmodel16$results$p_lt<=0.1,])
nrow(CHcoSmodel16$results[CHcoSmodel16$results$padj<=0.1&CHcoSmodel16$results$p_gt<=0.1,])
nrow(CHcoSmodel16$results[CHcoSmodel16$results$padj<=0.1&CHcoSmodel16$results$p_lt<=0.1,])
nrow(BWcoHmodel16$results[BWcoHmodel16$results$padj<=0.1&BWcoHmodel16$results$p_gt<=0.1,])
nrow(BWcoHmodel16$results[BWcoHmodel16$results$padj<=0.1&BWcoHmodel16$results$p_lt<=0.1,])
nrow(BWcoSmodel16$results[BWcoSmodel16$results$padj<=0.1&BWcoSmodel16$results$p_gt<=0.1,])
nrow(BWcoSmodel16$results[BWcoSmodel16$results$padj<=0.1&BWcoSmodel16$results$p_lt<=0.1,])
nrow(SPcoHmodel16$results[SPcoHmodel16$results$padj<=0.1&SPcoHmodel16$results$p_gt<=0.1,])
nrow(SPcoHmodel16$results[SPcoHmodel16$results$padj<=0.1&SPcoHmodel16$results$p_lt<=0.1,])
nrow(SPcoSmodel16$results[SPcoSmodel16$results$padj<=0.1&SPcoSmodel16$results$p_gt<=0.1,])
nrow(SPcoSmodel16$results[SPcoSmodel16$results$padj<=0.1&SPcoSmodel16$results$p_lt<=0.1,])


head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
head(CHcoHmodel$results[order(CHcoHmodel$results$p_gt,CHcoHmodel$results$padj),])
head(CHcoSmodel$results[order(CHcoSmodel$results$p_lt,CHcoSmodel$results$padj),])
head(CHcoSmodel$results[order(CHcoSmodel$results$p_gt,CHcoSmodel$results$padj),])
head(SPcoHmodel$results[order(SPcoHmodel$results$p_lt,SPcoHmodel$results$padj),])
head(SPcoHmodel$results[order(SPcoHmodel$results$p_gt,SPcoHmodel$results$padj),])
head(SPcoSmodel$results[order(SPcoSmodel$results$p_lt,SPcoSmodel$results$padj),])
head(SPcoSmodel$results[order(SPcoSmodel$results$p_gt,SPcoSmodel$results$padj),])
head(BWcoHmodel$results[order(BWcoHmodel$results$p_lt,BWcoHmodel$results$padj),])
head(BWcoHmodel$results[order(BWcoHmodel$results$p_gt,BWcoHmodel$results$padj),])
head(BWcoSmodel$results[order(BWcoSmodel$results$p_lt,BWcoSmodel$results$padj),])
head(BWcoSmodel$results[order(BWcoSmodel$results$p_gt,BWcoSmodel$results$padj),])

X <- head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_lt,CHcoHmodel16$results$padj),])
X <- rbind(X,head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_gt,CHcoHmodel16$results$padj),]))
X <- rbind(X,head(CHcoSmodel16$results[order(CHcoSmodel16$results$p_lt,CHcoSmodel16$results$padj),]))
X <- rbind(X,head(CHcoSmodel16$results[order(CHcoSmodel16$results$p_gt,CHcoSmodel16$results$padj),]))
X <- rbind(X,head(SPcoHmodel16$results[order(SPcoHmodel16$results$p_lt,SPcoHmodel16$results$padj),]))
X <- rbind(X,head(SPcoHmodel16$results[order(SPcoHmodel16$results$p_gt,SPcoHmodel16$results$padj),]))
X <- rbind(X,head(SPcoSmodel16$results[order(SPcoSmodel16$results$p_lt,SPcoSmodel16$results$padj),]))
X <- rbind(X,head(SPcoSmodel16$results[order(SPcoSmodel16$results$p_gt,SPcoSmodel16$results$padj),]))
X <- rbind(X,head(BWcoHmodel16$results[order(BWcoHmodel16$results$p_lt,BWcoHmodel16$results$padj),]))
X <- rbind(X,head(BWcoHmodel16$results[order(BWcoHmodel16$results$p_gt,BWcoHmodel16$results$padj),]))
X <- rbind(X,head(BWcoSmodel16$results[order(BWcoSmodel16$results$p_lt,BWcoSmodel16$results$padj),]))
X <- rbind(X,head(BWcoSmodel16$results[order(BWcoSmodel16$results$p_gt,BWcoSmodel16$results$padj),]))


Y <- head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
Y <- rbind(Y,head(CHcoHmodel$results[order(CHcoHmodel$results$p_gt,CHcoHmodel$results$padj),]))
Y <- rbind(Y,head(CHcoSmodel$results[order(CHcoSmodel$results$p_lt,CHcoSmodel$results$padj),]))
Y <- rbind(Y,head(CHcoSmodel$results[order(CHcoSmodel$results$p_gt,CHcoSmodel$results$padj),]))
Y <- rbind(Y,head(SPcoHmodel$results[order(SPcoHmodel$results$p_lt,SPcoHmodel$results$padj),]))
Y <- rbind(Y,head(SPcoHmodel$results[order(SPcoHmodel$results$p_gt,SPcoHmodel$results$padj),]))
Y <- rbind(Y,head(SPcoSmodel$results[order(SPcoSmodel$results$p_lt,SPcoSmodel$results$padj),]))
Y <- rbind(Y,head(SPcoSmodel$results[order(SPcoSmodel$results$p_gt,SPcoSmodel$results$padj),]))
Y <- rbind(Y,head(BWcoHmodel$results[order(BWcoHmodel$results$p_lt,BWcoHmodel$results$padj),]))
Y <- rbind(Y,head(BWcoHmodel$results[order(BWcoHmodel$results$p_gt,BWcoHmodel$results$padj),]))
Y <- rbind(Y,head(BWcoSmodel$results[order(BWcoSmodel$results$p_lt,BWcoSmodel$results$padj),]))
Y <- rbind(Y,head(BWcoSmodel$results[order(BWcoSmodel$results$p_gt,BWcoSmodel$results$padj),]))



#===============================================================================
#       Speculation site
#===============================================================================

# PCA

filterFun=function(o,f){
	o<-prune_samples(sample_data(o)[[2]]=="Speculation",o)
	prune_taxa(rowSums(otu_table(o))>5,o)
}

mypca <- plotPCA(biomITS,design="1",returnData=T,filterFun=filterFun,calcFactors=geoMeans)
myfiltbiom <- prune_samples(sample_data(biomITS)[[2]]=="Speculation",biomITS)
myfiltbiom <-  prune_taxa(rowSums(otu_table(myfiltbiom ))>5,myfiltbiom )
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

pdf("spec_its.pdf")
plotOrd(df,sample_data(myfiltbiom),design="condition",dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
dev.off

mypca <- plotPCA(biom,design="1",returnData=T,filterFun=filterFun,calcFactors=geoMeans)
myfiltbiom <- prune_samples(sample_data(biom)[[2]]=="Speculation",biom)
myfiltbiom <-  prune_taxa(rowSums(otu_table(myfiltbiom ))>5,myfiltbiom )
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

pdf("spec_its.pdf")
plotOrd(df,sample_data(myfiltbiom),design="condition",dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
dev.off()


# differential analysis

site="Speculation"
dds <- phylo_to_des(biomITS ,fit=F,calcFactors=geoMeans)
dds <- dds[ rowSums(counts(dds)) > 5, ]
#dds <- dds[,dds$site==site]
dds$condition <- droplevels(dds$condition)
dds$site <- droplevels(dds$site)
dds$tree <- droplevels(dds$tree)

design=~condition#~tree + condition
design(dds) <- design
dds <- DESeq(dds,parallel=T)
alpha <- 0.05

contrast=c("condition","Healthy","Symptom")
res <- results(dds,contrast=contrast,alpha=alpha,parallel=T)

myfiltbiom <- prune_samples(sample_data(biomITS)$site==site,mybiomITS)

res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)

write.table(res.merge,"speculation.fun.res",sep="\t",quote=F,na="",row.names=F)


dds <- phylo_to_des(biom ,fit=F)
dds <- dds[ rowSums(counts(dds)) > 5, ]
#dds <- dds[,dds$site==site]
dds$condition <- droplevels(dds$condition)
dds$site <- droplevels(dds$site)
dds$tree <- droplevels(dds$tree)

design=~condition#~tree + condition
design(dds) <- design
dds <- DESeq(dds,parallel=T)
alpha <- 0.05

contrast=c("condition","Healthy","Symptom")
res <- results(dds,contrast=contrast,alpha=alpha,parallel=T)

myfiltbiom <- prune_samples(sample_data(biom)$site==site,biom)

res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)

write.table(res.merge,"speculation.bac.res",sep="\t",quote=F,na="",row.names=T)


## Volcano plots

pdf("volcano_spec.pdf")

with(res.merge.its,plot(log2FoldChange,log10(baseMean),pch=20, xlim=c(-2.5,2)))
with(subset(res.merge.its, padj<0.05 ), points(log2FoldChange, log10(baseMean), pch=20, col="red"))
with(subset(res.merge.its, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="orange"))
with(subset(res.merge.its, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="green"))


dev.off()






### Other stuff


# using combined biological replicates (I think this is a bad idea..)
groupby=paste(dds$condition,dds$tree,sep="_")
dds2 <- collapseReplicates2(dds,groupby=groupby) # uses mean for grouping rather than sum
design(dds2) <- ~condition
dds2 <- DESeq(dds2,parallel=T)
res2 <- results(dds2,alpha=alpha,parallel=T)

pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$tree))

df <- t(data.frame(t(mypca$x)*mypca$percentVar))
d <- t(data.frame(t(pc.res)*mypca$percentVar))

sample_data(myfiltbiom)$loc <- as.numeric(levels(sample_data(myfiltbiom)$tree)[sample_data(myfiltbiom)$tree])
plotOrd(df,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=1,dimy=2,xlabel="PC1",ylabel="PC2")
plotOrd(d,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=1,dimy=2,xlabel="PC1",ylabel="PC2")
plotOrd(df,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=2,dimy=3,xlabel="PC2",ylabel="PC3")
plotOrd(d,sample_data(myfiltbiom),shapes="condition",design="loc",continuous=T,dimx=2,dimy=3,xlabel="PC2",ylabel="PC3")


pc.res <- resid(aov(mypca$x~sample_data(myfiltbiom)$tree))

mybiom <- biom
#mybiom <- biomITS

myfiltbiom <-  prune_taxa(rowSums(otu_table(mybiom ))>5,mybiom )
mypca <- plotPCA(myfiltbiom ,design="1",returnData=T,calcFactors=geoMeans,filterFun=filterFun,filter=...)


#mypca <- plotPCA(myfiltbiom ,design="1",returnData=T,calcFactors=geoMeans) # ITS only (bypass log geometric means error)


df <- t(data.frame(t(mypca$x)*mypca$percentVar))
pdf("all_S_pca.pdf")
#pdf("all_ITS_pca.pdf")

plotOrd(df,sample_data(myfiltbiom),design=c("site","condition"),dimx=1,dimy=2,xlabel="PC1",ylabel="PC2",pointSize=1.5,cbPallete=T)
plotOrd(df,sample_data(myfiltbiom),design=c("site","condition"),dimx=2,dimy=3,xlabel="PC2",ylabel="PC3",pointSize=1.5,cbPallete=T)
dev.off()
