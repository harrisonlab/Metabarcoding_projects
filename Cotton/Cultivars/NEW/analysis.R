#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(2))
library(data.table)
library(dplyr)
library(plyr)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
library(vegan)
library(ape)
library(GUniFrac)
library(ggplot2)
library(reshape2)
library(tibble)
library(lmPerm)


#===============================================================================
#       Load data
#===============================================================================

ubiom_BAC_Rhiz <- loadData("final_16s_Rhiz_otu_table.txt",
													 "16S cultivar meta data.csv",
													 "16S_Rhiz.taxa",
													 "16S_Rhiz.phy",
													 RHB="BAC_RHIZ")
ubiom_BAC_Endo <- loadData("final_16s_Endo_otu_table.txt",
													 "16S cultivar meta data.csv",
													 "16S_Endo.taxa",
													 RHB="BAC_ENDO")
ubiom_FUN_Rhiz      <- loadData("FUN_RHIZ.otus_table.txt",
													 "ITS cultivar meta data.csv",
													 "FUN_RHIZ.taxa",
													 "FUN_RHIZ.phy",
													 RHB="FUN_RHIZ")
ubiom_FUN_Endo      <- loadData("FUN_ENDO.otus_table.txt",
													 "ITS cultivar meta data.csv",
													 "FUN_ENDO.taxa",
													 "FUN_ENDO.phy",
													 RHB="FUN_ENDO")

colnames(ubiom_FUN_Rhiz$countData) <- gsub("(DMP[0-9]*\\.)([A-Z]*[0-9]*)(\\.)([ER][0-9])(_.*)","\\2\\4", colnames(ubiom_FUN_Rhiz$countData))
colnames(ubiom_FUN_Endo$countData) <- gsub("(DMP[0-9]*\\.)([A-Z]*[0-9]*)(\\.)([ER][0-9])(_.*)","\\2\\4", colnames(ubiom_FUN_Endo$countData))

ubiom_FUN_Rhiz$countData <- ubiom_FUN_Rhiz$countData[,grep("R",colnames(ubiom_FUN_Rhiz$countData))]
ubiom_FUN_Rhiz$colData	 <- ubiom_FUN_Rhiz$colData[grep("R",rownames(ubiom_FUN_Rhiz$colData)),2:3, drop = FALSE]
ubiom_FUN_Endo$countData <- ubiom_FUN_Endo$countData[,grep("E",colnames(ubiom_FUN_Endo$countData))]
ubiom_FUN_Endo$colData	 <- ubiom_FUN_Endo$colData[grep("E",rownames(ubiom_FUN_Endo$colData)),2:3, drop = FALSE]
ubiom_FUN_Endo$countData <- ubiom_FUN_Endo$countData[rownames(ubiom_FUN_Endo$countData)!="OTU1",]
ubiom_FUN_Endo$taxData <- ubiom_FUN_Endo$taxData[rownames(ubiom_FUN_Endo$taxData)!="OTU1",]
ubiom_FUN_Rhiz$countData <- ubiom_FUN_Rhiz$countData[,c(-1,-3,-27)]

colnames(ubiom_FUN_Rhiz$countData) <- gsub("\\..*","", colnames(ubiom_FUN_Rhiz$countData))


# and the metadata/countData naming convention in BAC rhiz and endo is different
rownames(ubiom_BAC_Rhiz$colData) <- sub("E","R",rownames(ubiom_BAC_Rhiz$colData))

#===============================================================================
#       Create DEseq objects
#===============================================================================

ubiom_BAC_Rhiz$dds <- ubiom_to_des(ubiom_BAC_Rhiz)
ubiom_BAC_Endo$dds <- ubiom_to_des(ubiom_BAC_Endo)
ubiom_FUN_Rhiz$dds <- ubiom_to_des(ubiom_FUN_Rhiz)
ubiom_FUN_Endo$dds <- ubiom_to_des(ubiom_FUN_Endo)

ubiom_all <- list(ubiom_BAC_Rhiz=ubiom_BAC_Rhiz,ubiom_BAC_Endo=ubiom_BAC_Endo,ubiom_FUN_Rhiz=ubiom_FUN_Rhiz, ubiom_FUN_Endo= ubiom_FUN_Endo)
ubiom_all_bak <- ubiom_all

#===============================================================================
#       filter data
#===============================================================================

# filter samples for less than 1000 reads and suseptibility
samfilter <- lapply(ubiom_all,function(u) colSums(counts(u$dds))>=1000&(colData(u$dds)$Susceptibility>35|colData(u$dds)$Susceptibility<10))
# filter OTUs
otufilter <- lapply(ubiom_all,function(u) rowSums(counts(u$dds, normalize=T))>4)

ubiom_BAC_Rhiz$dds <- ubiom_BAC_Rhiz$dds[otufilter$ubiom_BAC_Rhiz,samfilter$ubiom_BAC_Rhiz]
ubiom_BAC_Endo$dds <- ubiom_BAC_Endo$dds[otufilter$ubiom_BAC_Endo,samfilter$ubiom_BAC_Endo]
ubiom_FUN_Rhiz$dds <- ubiom_FUN_Rhiz$dds[otufilter$ubiom_FUN_Rhiz,samfilter$ubiom_FUN_Rhiz]
ubiom_FUN_Endo$dds <- ubiom_FUN_Endo$dds[otufilter$ubiom_FUN_Endo,samfilter$ubiom_FUN_Endo]

# filter taxData
ubiom_BAC_Rhiz$taxData <- ubiom_BAC_Rhiz$taxData[which(rownames(ubiom_BAC_Rhiz$taxData)%in%rownames(counts(ubiom_BAC_Rhiz$dds))),]
ubiom_BAC_Endo$taxData <- ubiom_BAC_Endo$taxData[which(rownames(ubiom_BAC_Endo$taxData)%in%rownames(counts(ubiom_BAC_Endo$dds))),]
ubiom_FUN_Rhiz$taxData <- ubiom_FUN_Rhiz$taxData[which(rownames(ubiom_FUN_Rhiz$taxData)%in%rownames(counts(ubiom_FUN_Rhiz$dds))),]
ubiom_FUN_Endo$taxData <- ubiom_FUN_Endo$taxData[which(rownames(ubiom_FUN_Endo$taxData)%in%rownames(counts(ubiom_FUN_Endo$dds))),]

ubiom_FUN_Endo$dds <- ubiom_FUN_Endo$dds[rownames(ubiom_FUN_Endo$dds)!="OTU1",]
ubiom_FUN_Endo$taxData <- ubiom_FUN_Endo$taxData[rownames(ubiom_FUN_Endo$taxData)!="OTU1",]

#===============================================================================
#       Plot/table frequencies
#===============================================================================

# bacteria
invisible(mapply(assign, names(ubiom_BAC_Rhiz), ubiom_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
rhiz <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.8,proportional=T)
invisible(mapply(assign, names(ubiom_BAC_Endo), ubiom_BAC_Endo, MoreArgs=list(envir = globalenv())))
endo <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.8,proportional=T)
md1 <- melt(rhiz,id=colnames(rhiz)[1])
md2 <- melt(endo,id=colnames(endo)[1])
levels(md1$variable)[1] <- "Rhizosphere"
levels(md2$variable)[1] <- "Endophyte"

md1$phylum  <- factor(md1$phylum, levels=unique(md1$phylum[order(md1$value,decreasing=T)]),ordered=T)
md2$phylum  <- factor(md2$phylum, levels=unique(md2$phylum[order(md2$value,decreasing=T)]),ordered=T)

md2 <- rbind(md1,md2)
#md2$variable <- factor(md2$phylum, levels=levels(md2$phylum)[order(levels(md2$value))]  )
md2$value <- as.numeric(md2$value)
md2$Samples <- md2$variable

g <- ggplot(md2[md2$value>0.3,],aes(x=phylum,y=value,fill=Samples)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
#g <- g + geom_bar(stat="identity",colour="white",position = "dodge")
g <- g + scale_fill_manual(values = c("Rhizosphere" = "black", "Endophyte" = "orange"))
g <- g  + xlab("")+ ylab("Frequency (%)")
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
	axis.title.y=element_text(size=(14-2)))
ggsave("BAC_frequency.pdf",g,width=8,height=7)

invisible(mapply(assign, names(ubiom_BAC_Rhiz), ubiom_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
x1 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.65,proportional=T,taxon="phylum")
invisible(mapply(assign, names(ubiom_BAC_Endo), ubiom_BAC_Endo, MoreArgs=list(envir = globalenv())))
x2 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.65,proportional=T,taxon="phylum")
x1$all <- as.numeric(x1$all)
x2$all <- as.numeric(x2$all)
colnames(x1)[2] <- "Rhizosphere"
colnames(x2)[2] <- "Endophyte"
write.table(full_join(x1,x2),"BAC_phylum_freq.csv",sep=",",row.names=F,quote=F,na="")

invisible(mapply(assign, names(ubiom_BAC_Rhiz), ubiom_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
x1 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.65,proportional=T,taxon="class")
invisible(mapply(assign, names(ubiom_BAC_Endo), ubiom_BAC_Endo, MoreArgs=list(envir = globalenv())))
x2 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.65,proportional=T,taxon="class")
x1$all <- as.numeric(x1$all)
x2$all <- as.numeric(x2$all)
colnames(x1)[2] <- "Rhizosphere"
colnames(x2)[2] <- "Endophyte"
write.table(full_join(x1,x2),"BAC_class_freq.csv",sep=",",row.names=F,quote=F,na="")

invisible(mapply(assign, names(ubiom_BAC_Rhiz), ubiom_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
x1 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="order")
invisible(mapply(assign, names(ubiom_BAC_Endo), ubiom_BAC_Endo, MoreArgs=list(envir = globalenv())))
x2 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="order")
x1$all <- as.numeric(x1$all)
x2$all <- as.numeric(x2$all)
colnames(x1)[2] <- "Rhizosphere"
colnames(x2)[2] <- "Endophyte"
write.table(full_join(x1,x2),"BAC_order_freq.csv",sep=",",row.names=F,quote=F,na="")

# fungi

invisible(mapply(assign, names(ubiom_FUN_Rhiz), ubiom_FUN_Rhiz, MoreArgs=list(envir = globalenv())))
rhiz <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.65,proportional=T,taxon="class")
invisible(mapply(assign, names(ubiom_FUN_Endo), ubiom_FUN_Endo, MoreArgs=list(envir = globalenv())))
endo <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.65,proportional=T,taxon="class")
md1 <- melt(rhiz,id=colnames(rhiz)[1])
md2 <- melt(endo,id=colnames(endo)[1])
levels(md1$variable)[1] <- "Rhizosphere"
levels(md2$variable)[1] <- "Endophyte"

#keep <- md1$value>=0.3|md2$value>=0.3

md1$class  <- factor(md1$class, levels=unique(md1$class[order(md1$value,decreasing=T)]),ordered=T)
md2$class  <- factor(md2$class, levels=unique(md2$class[order(md2$value,decreasing=T)]),ordered=T)

levels(md1$class) <- gsub("_cls.*"," (Ins sed)",levels(md1$class))
levels(md2$class) <- gsub("_cls.*"," (Ins sed)",levels(md2$class))

md2 <- rbind(md1,md2)
#md2$class <- factor(md2$class, levels=levels(md2$class)[order(md2$value,decreasing=F)])
md2$value <- as.numeric(md2$value)
md2$Samples <- md2$variable

g <- ggplot(md2[md2$value>0.3,],aes(x=class,y=value,fill=Samples)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
#g <- g + geom_bar(stat="identity",colour="white",position = "dodge")
g <- g + scale_fill_manual(values = c("Rhizosphere" = "black", "Endophyte" = "orange"))
g <- g  + xlab("")+ ylab("Frequency (%)")
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
	axis.title.y=element_text(size=(14-2)))
ggsave("FUN_frequency.pdf",g,width=10,height=7)


invisible(mapply(assign, names(ubiom_FUN_Rhiz), ubiom_FUN_Rhiz, MoreArgs=list(envir = globalenv())))
x1 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="class")
invisible(mapply(assign, names(ubiom_FUN_Endo), ubiom_FUN_Endo, MoreArgs=list(envir = globalenv())))
x2 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="class")
x1$all <- as.numeric(x1$all)
x2$all <- as.numeric(x2$all)
colnames(x1)[2] <- "Rhizosphere"
colnames(x2)[2] <- "Endophyte"
write.table(full_join(x1,x2),"FUN_class_freq.csv",sep=",",row.names=F,quote=F)

invisible(mapply(assign, names(ubiom_FUN_Rhiz), ubiom_FUN_Rhiz, MoreArgs=list(envir = globalenv())))
x1 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="order")
invisible(mapply(assign, names(ubiom_FUN_Endo), ubiom_FUN_Endo, MoreArgs=list(envir = globalenv())))
x2 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="order")
x1$all <- as.numeric(x1$all)
x2$all <- as.numeric(x2$all)
colnames(x1)[2] <- "Rhizosphere"
colnames(x2)[2] <- "Endophyte"
write.table(full_join(x1,x2),"FUN_order_freq.csv",sep=",",row.names=F,quote=F)


invisible(mapply(assign, names(ubiom_FUN_Rhiz), ubiom_FUN_Rhiz, MoreArgs=list(envir = globalenv())))
x1 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="family")
invisible(mapply(assign, names(ubiom_FUN_Endo), ubiom_FUN_Endo, MoreArgs=list(envir = globalenv())))
x2 <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.5,proportional=T,taxon="family")
x1$all <- as.numeric(x1$all)
x2$all <- as.numeric(x2$all)
colnames(x1)[2] <- "Rhizosphere"
colnames(x2)[2] <- "Endophyte"
write.table(full_join(x1,x2),"FUN_family_freq.csv",sep=",",row.names=F,quote=F)

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects - run only one
invisible(mapply(assign, names(ubiom_BAC_Rhiz), ubiom_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC_Endo), ubiom_BAC_Endo, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_FUN_Rhiz), ubiom_FUN_Rhiz, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_FUN_Endo), ubiom_FUN_Endo, MoreArgs=list(envir = globalenv())))

# add various susceptibility bins
dds$Sus <- "Susceptible"
dds$Sus[dds$Susceptibility<10] <- "Tolerant"
#dds$Sus[dds$Susceptibility>10&dds$Susceptibility<30] <- "Semi-tolerant"
dds$Sus <- as.factor(dds$Sus)

#===============================================================================
#       Alpha diversity analysis - RUN BEFORE FILTERING OUT ANY LOW COUNT OTUS
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),colour="Cultivar",design="Susceptibility",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(0,5000,"Chao1")))

### permutation based anova on diversity index ranks ###

# get alpha diversity indices
alpha_ord <- plot_alpha(counts(dds,normalize=T),colData(dds),design="Susceptibility",colour="Cultivar",returnData=T)

# join diversity indices and metadata
alpha_ord <- 	as.data.table(left_join(alpha_ord,as.data.frame(colData(dds))%>% mutate(Samples = rownames(colData(dds)))))

# perform anova for each index  - need to work out what we need to test for

sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
	setkey(alpha_ord,S.chao1)
	print("Chao1")
	summary(aovp(as.numeric(as.factor(alpha_ord$S.chao1))~Sus,alpha_ord))
	setkey(alpha_ord,shannon)
	print("Shannon")
	summary(aovp(as.numeric(as.factor(alpha_ord$shannon))~Sus,alpha_ord))
	setkey(alpha_ord,simpson)
	print("simpson")
	summary(aovp(as.numeric(as.factor(alpha_ord$simpson))~Sus,alpha_ord))
sink()


#===============================================================================
#       Beta diversity
#===============================================================================

### NMDS ###

# phyloseq has functions (using Vegan) for making NMDS plots
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,as.data.frame(colData(dds))))

# add tree to phyloseq object
phy_tree(myphylo) <- nj(as.dist(phylipData))

# calculate NMDS ordination using weighted unifrac scores
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

# plot with plotOrd (or use plot_ordination)
ggsave(paste(RHB,"Unifrac_NMDS.pdf",sep="_"),plotOrd(ordu$points,colData(dds),design="Cultivar",shape="Sus",xlabel="NMDS1",ylabel="NMDS2",alpha=0.75),width=7,height=7)

# permanova of unifrac distance
sink(paste(RHB,"PERMANOVA_unifrac.txt",sep="_"))
 print("weighted")
 adonis(distance(myphylo,"unifrac",weighted=T)~Sus,colData(dds),parallel=12,permutations=9999)
 print("unweighted")
 adonis(distance(myphylo,"unifrac",weighted=F)~Sus,colData(dds),parallel=12,permutations=9999)
sink()

#===============================================================================
#       PCA
#===============================================================================

# perform PC decompossion on DES object
mypca <- des_to_pca(dds)

### ANOVA ###

# write to file
sink(paste(RHB,"PCA.txt"))

cat("
# Variance in first 4 PC scores \n")
round(mypca$percentVar[1:4] * 100,2)

cat("
# ANOVA of first four PC scores \n")
apply(mypca$x[,1:4],2,function(x) summary(aov(x~Cultivar,as.data.frame(dds@colData))))

# get sum of squares for all PC scores
sum_squares <- 

# name sum_squares columns
colnames(sum_squares) <- c("Cultivar","residual")

# proportion of total sum of squares for PC
perVar <- t(apply(sum_squares,1,prop.table)) * mypca$percentVar

# check - should equal 1
# sum(colSums(perVar))

cat("t(apply(mypca$x,2,function(x)t(summary(aov(x~Cultivar,as.data.frame(dds@colData)))[[1]][2]))) - 
  t(apply(mypca$x,2,function(x)t(summary(aov(x~Sus,as.data.frame(dds@colData)))[[1]][2])))
# % total variance explained by the aov model and residual \n")
colSums(perVar)/sum(colSums(perVar))*100

# end write to file
sink()

### Plot ###

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
DF  <- t(data.frame(t(mypca$x)*mypca$percentVar))

# output pdf
pdf(paste(RHB,"PCA.pdf",sep="_"))

# plot PC1 vs PC2
plotOrd(DF,dds@colData,design="Cultivar",shape="Sus",axes=c(1,2),alpha=0.75)

# plot PC2 vs PC3
plotOrd(DF,dds@colData,design="Cultivar",shape="Sus",axes=c(2,3),alpha=0.75)

# plot PC3 vs PC4
plotOrd(DF,dds@colData,design="Cultivar",shape="Sus",axes=c(3,4),alpha=0.75)

# write to file
dev.off()

# output pdf
pdf(paste(RHB,"PCA_susceptibility.pdf",sep="_"))

# plot PC1 vs PC2
plotOrd(DF,dds@colData,design="Susceptibility",axes=c(1,2),alpha=0.75,continuous=T)

# plot PC2 vs PC3
plotOrd(DF,dds@colData,design="Susceptibility",axes=c(2,3),alpha=0.75,continuous=T)

# plot PC3 vs PC4
plotOrd(DF,dds@colData,design="Susceptibility",axes=c(3,4),alpha=0.75,continuous=T)

# write to file
dev.off()

#===============================================================================
#       RDA plots
#===============================================================================

# agregate counts at phylum level
combinedTaxa <- combineTaxa2(taxData[,-8],rank="phylum",confidence=0.8,returnFull=T)
phylumData <- combCounts(combinedTaxa,counts(dds,normalize=T))
phylumTaxa <- combTaxa(combinedTaxa,taxData[,-8])
phylumData <- aggregate(phylumData,by=list(taxaConfVec(phylumTaxa,conf=0.8,level=2)),sum)
rownames(phylumData) <- phylumData[,1]
phylumData <- phylumData[,-1]

myrda <- rda(t(phylumData)~Cultivar,colData(dds))

species <- as.data.frame(scores(myrda,scaling="symmetric")$species)
species$phylum<-rownames(species)
# filter out phyla with less than 1% of total abundance
species <- species[((rowSums(phylumData)/sum(phylumData))*100)>=2,]
#centroids <- centroids[centroids$phylum %in% rhiz$phylum[rhiz$all>=4],]
sites   <- scores(myrda)$sites

ggsave(paste(RHB,"RDA_phylum_biplot.pdf",sep="_"),plotOrd(sites,dds@colData,design="Cultivar",shape="Sus",axes=c(1,2),alpha=0.75,legend=T,continuous=F,ylims=c(-100,100),xlims=c(-100,100))  +
 # geom_point(data=species,aes(x=RDA1,y=RDA2,fill=phylum,shape=phylum), size=4,inherit.aes=F) +
	geom_segment(data=species,aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
	geom_label(data=species,aes(label=phylum,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))

ggsave(paste(RHB,"RDA_phylum_biplot_continuous.pdf",sep="_"),plotOrd(sites,dds@colData,design="Susceptibility",axes=c(1,2),alpha=0.75,legend=T,continuous=T,ylims=c(-100,100),xlims=c(-100,100))  +
 # geom_point(data=species,aes(x=RDA1,y=RDA2,fill=phylum,shape=phylum), size=4,inherit.aes=F) +
	geom_segment(data=species,aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
	geom_label(data=species,aes(label=phylum,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))


# agregate counts at class level
combinedTaxa <- combineTaxa2(taxData[,-8],rank="class",confidence=0.8,returnFull=T)
classData <- combCounts(combinedTaxa,counts(dds,normalize=T))
classTaxa <- combTaxa(combinedTaxa,taxData[,-8])
classData <- aggregate(classData,by=list(taxaConfVec(classTaxa,conf=0.8,level=3)),sum)
rownames(classData) <- classData[,1]
classData <- classData[,-1]

myrda <- rda(t(classData)~Cultivar,colData(dds))

species <- as.data.frame(scores(myrda,scaling="symmetric")$species)
species$class<-rownames(species)
# filter out phyla with less than 1% of total abundance
species <- species[((rowSums(classData)/sum(classData))*100)>=5,]
#centroids <- centroids[centroids$phylum %in% rhiz$phylum[rhiz$all>=4],]
sites   <- scores(myrda)$sites

ggsave(paste(RHB,"RDA_class_biplot.pdf",sep="_"),plotOrd(sites,dds@colData,design="Cultivar",shape="Sus",axes=c(1,2),alpha=0.75,legend=T,ylims=c(-100,100),xlims=c(-100,100))+#,ylims=c(-140,140))  +#,xlims=c(-150,150)
 # geom_point(data=species,aes(x=RDA1,y=RDA2,fill=phylum,shape=phylum), size=4,inherit.aes=F) +
	geom_segment(data=species,aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
	geom_label(data=species,aes(label=class,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))

ggsave(paste(RHB,"RDA_class_biplot_continuous.pdf",sep="_"),plotOrd(sites,dds@colData,design="Susceptibility",axes=c(1,2),alpha=0.75,legend=T,continuous=T)+#,ylims=c(-140,140))  +#,xlims=c(-150,150)
 # geom_point(data=species,aes(x=RDA1,y=RDA2,fill=phylum,shape=phylum), size=4,inherit.aes=F) +
	geom_segment(data=species,aes(x=0,y=0,xend=RDA1,yend=RDA2), size=0.5,arrow=arrow(),inherit.aes=F) +
	geom_label(data=species,aes(label=class,x=(RDA1/2),y=(RDA2/2)),size=2,inherit.aes=F))

# filter fungi
# rhiz 2 and 3
# endo 2 and 1
# bac rhiz 2 and 4
# bac endo 1 and 1

#===============================================================================
#       Statistical analysis
#===============================================================================


# p value for FDR cutoff
alpha <- 0.1

# drop any unused levels in Cultivar
dds@colData <- droplevels(dds@colData)

# add model to the DES object
design(dds) <- ~Sus

# calculate fit
dds <- DESeq(dds,parallel=T)

# Tolerant vs suspetible
res <- results(dds,alhpa=alpha,parallel=T,contrast=c("Sus","Tolerant","Susceptible"))

# Tolerant vs Semi-tolernat
#res2 <- results(dds,alhpa=alpha,parallel=T,contrast=c("Sus","Tolerant","Semi-tolerant"))

# Semi-tolernat vs suspetible
#res3 <- results(dds,alhpa=alpha,parallel=T,contrast=c("Sus","Semi-tolerant","Susceptible"))

# Tolerant vs Semi-tolernat + suspetible
#res4 <- results(dds,alhpa=alpha,parallel=T,contrast=c(1,0,1))

#dds2 <- dds
#levels(dds2$Sus)[levels(dds2$Sus)=="Semi-tolerant"] <- "Susceptible"
#design(dds2) <- ~Sus
#dds2 <- DESeq(dds2,parallel=T)
#res4 <- results(dds2,alhpa=alpha,parallel=T)

fwrite(data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))),paste0(RHB,"_Tolerant_vs_Suspetible.txt"),sep="\t",quote=F,na="")
#fwrite(data.table(inner_join(data.table(OTU=rownames(res2),as.data.frame(res2)),data.table(OTU=rownames(taxData),taxData))),paste0(RHB,"_Tolerant_vs_Semi-tolerant.txt"),sep="\t",quote=F,na="")
#fwrite(data.table(inner_join(data.table(OTU=rownames(res3),as.data.frame(res3)),data.table(OTU=rownames(taxData),taxData))),paste0(RHB,"_Semi-tolerant_vs_Suspetible.txt"),sep="\t",quote=F,na="")
#fwrite(data.table(inner_join(data.table(OTU=rownames(res4),as.data.frame(res4)),data.table(OTU=rownames(taxData),taxData))),paste0(RHB,"_Tolerant_vs_rest.txt"),sep="\t",quote=F,na="")

## MA plot

ggsave(paste(RHB,"ma_plot.pdf",sep="_"),plot_ma(res[,c(2,1,5)],crush=T,legend=T))
