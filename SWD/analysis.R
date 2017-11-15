#===============================================================================
#       Load libraries
#===============================================================================

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")

setEPS() # sets postscript defaults

#===============================================================================
#       Load data 
#===============================================================================

mybiom <- import_biom("first/16S.taxa.biom") 
sample_data(mybiom) <- read.table("first/colData",header=T,sep="\t",row.names=1)
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65,level=7)
sample_data(mybiom)$all <- "all"
mybiom_fb <- mybiom


mybiom <- import_biom("second/16S.taxa.biom")
sample_data(mybiom) <- read.table("second/colData",header=T,sep="\t",row.names=1)
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65,7)
sample_data(mybiom)$all <- "all"
mybiom_sb <- mybiom

 
mybiom <- merge_phyloseq(import_biom("first/ITS1.taxa.biom"),import_biom("first/ITS2.taxa.biom"))
sample_data(mybiom) <- read.table("first/colData",header=T,sep="\t",row.names=1)
tax_table(mybiom)[,8] <- 1
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65,6)
sample_data(mybiom)$all <- "all"
mybiom_ff <- mybiom

mybiom <- merge_phyloseq(import_biom("second/ITS1.taxa.biom"),import_biom("second/ITS2.taxa.biom"))
sample_data(mybiom) <- read.table("second/colData",header=T,sep="\t",row.names=1)
tax_table(mybiom)[,8] <- 1
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom),0.65,6)
sample_data(mybiom)$all <- "all"
mybiom_sf <- mybiom


#===============================================================================
#       diversity analysis
#===============================================================================


mybiom <- mybiom_sb
dds <- phylo_to_des(mybiom)

#mybiom@otu_table@.Data <- round(counts(dds,normalize=T),0)

# need to do this as round is flooring some of the reads down to zero - this will affect Chao1 calculations
mybiom@otu_table@.Data <- counts(dds,normalize=T)
mybiom@otu_table@.Data[mybiom@otu_table@.Data==0] <- NA
mybiom@otu_table@.Data[mybiom@otu_table@.Data<1] <- 1
mybiom@otu_table@.Data[is.na(mybiom@otu_table@.Data)] <- 0
mybiom@otu_table@.Data <- round(mybiom@otu_table@.Data,0)
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

postscript("alpha_colour_final_v4.eps",height=28,width=28)
plot_richness(mybiom,x="condition",color="Sex",measures=c("Chao1", "Shannon", "Simpson"),textSize=56,size=12,limits=c(0,700)) +labs(x="Sex")+scale_colour_manual(values=cbbPalette)
dev.off()


postscript("alpha_bw_final_v3.eps",height=12.5,width=12.5)
plot_richness(mybiom,x="condition",shape="Sex",measures=c("Chao1", "Shannon", "Simpson"),textSize=28,size=3.5,limits=c(0,700))+labs(x="Sex")
dev.off()


postscript("alpha_bw_Poster.eps",height=32,width=32)
plot_richness(mybiom,x="condition",shape="Sex",measures=c("Chao1", "Shannon", "Simpson"),textSize=28,size=3.5)+labs(x="Sex")
dev.off()



#ANOVA

mybiom <- mybiom_sb
dds <- phylo_to_des(mybiom)

#mybiom@otu_table@.Data <- round(counts(dds,normalize=T),0)

# need to do this as round is flooring some of the reads down to zero - this will affect Chao1 calculations
mybiom@otu_table@.Data <- counts(dds,normalize=T)
mybiom@otu_table@.Data[mybiom@otu_table@.Data==0] <- NA
mybiom@otu_table@.Data[mybiom@otu_table@.Data<1] <- 1
mybiom@otu_table@.Data[is.na(mybiom@otu_table@.Data)] <- 0
mybiom@otu_table@.Data <- round(mybiom@otu_table@.Data,0)
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)
mybiom <- prune_samples(sample_data(mybiom)$Sex!="FM",mybiom)


df <- plot_richness(mybiom,x="condition",shape="Sex",measures=c("Chao1", "Shannon", "Simpson"),returnData=T)
summary(aov(Chao1~Sex*condition,df))
summary(aov(Shannon~Sex*condition,df))
summary(aov(Simpson~Sex*condition,df))

options(contrasts = c("contr.sum","contr.poly"))
#options(contrasts = c("contr.treatment","contr.poly"))

model <- lm(Chao1~condition*Sex,df)
drop1(model, .~., test="F")


mybiom <- mybiom_sf
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <-counts(dds,normalize=T)
mybiom <- prune_samples(sample_data(mybiom)$Sex!="FM",mybiom)
mybiom@otu_table@.Data[mybiom@otu_table@.Data==0] <- NA
mybiom@otu_table@.Data[mybiom@otu_table@.Data<1] <- 1
mybiom@otu_table@.Data[is.na(mybiom@otu_table@.Data)] <- 0
mybiom@otu_table@.Data <- round(mybiom@otu_table@.Data,0)
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)
mybiom <- prune_samples(sample_data(mybiom)$Sex!="FM",mybiom)


df <- plot_richness(mybiom,x="condition",shape="Sex",measures=c("Chao1", "Shannon", "Simpson"),returnData=T)
summary(aov(Chao1~condition*Sex,df))
summary(aov(Shannon~condition*Sex,df))
summary(aov(Simpson~condition*Sex,df))


postscript("alpha_poster.eps",height=28,width=28)
plot_richness(mybiom,x="condition",color="Sex",measures=c("Chao1", "Shannon", "Simpson"),textSize=58,size=10) +labs(x="Sex")+scale_colour_manual(values=cbbPalette)
dev.off()

### Beta
mybiom <- mybiom_sb
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

adonis(t(as.data.frame(as.matrix(otu_table(myfiltbiom))))~condition*Sex,as.data.frame(as.matrix(sample_data(myfiltbiom))),parallel=12,permutations=9999)

mybiom <- mybiom_sf
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)
dds <- phylo_to_des(mybiom,fitType="local")
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

adonis(t(as.data.frame(as.matrix(otu_table(myfiltbiom))))~condition*Sex,as.data.frame(as.matrix(sample_data(myfiltbiom))),parallel=12,permutations=9999)


#===============================================================================
#       taxonomy plots
#===============================================================================

### bacteria

## 2015

mybiom <- mybiom_fb
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
tax_table(myfiltbiom)[,15] <- sub("\\(g\\)","",tax_table(myfiltbiom)[,15])
postscript("b_y1_all_final2.eps",width=10,height=10)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,bw=T,textSize=22,trans=F,legend=F,ylims=c(0,60),margins=c(0.5,0.2,0.2,1.5),ylab="% normalised reads")
dev.off()

## 2016

mybiom <- mybiom_sb
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
#dds <- phylo_to_des(myfiltbiom)
#vst <-  varianceStabilizingTransformation(dds)
#myfiltbiom@otu_table@.Data <- assay(vst)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)$condition=="N",myfiltbiom)
tax_table(myfiltbiom)[,15] <- sub("\\(g\\)","",tax_table(myfiltbiom)[,15])
tax_table(myfiltbiom) <- sub("_"," ",tax_table(myfiltbiom))
postscript("b_y2_N_final2.eps",width=10,height=10)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,bw=T,textSize=22,trans=F,legend=F,ylims=c(0,60),margins=c(0.5,0.2,0.2,1.5),ylab="% normalised reads")
dev.off()

### fungi

## 2015

mybiom <- mybiom_ff
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
tax_table(myfiltbiom)[,15] <- sub("\\(g\\)","",tax_table(myfiltbiom)[,15])
tax_table(myfiltbiom) <- sub("_fam_I.*\\(","(",tax_table(myfiltbiom))
postscript("f_y1_all_final2.eps",width=10,height=10)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,bw=T,textSize=22,trans=F,ylab="% normalised reads")#+theme(text=element_text(size=28))+scale_y_continuous(labels=scaleFUN)
dev.off()

## 2016

mybiom <- mybiom_sf
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
#dds <- phylo_to_des(myfiltbiom)
#vst <-  varianceStabilizingTransformation(dds)
#myfiltbiom@otu_table@.Data <- assay(vst)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)$condition=="N",myfiltbiom)
tax_table(myfiltbiom)[,15] <- sub("\\(g\\)","",tax_table(myfiltbiom)[,15])
tax_table(myfiltbiom) <- sub("_fam_I.*\\(","(",tax_table(myfiltbiom))
postscript("f_y2_N_final2.eps",width=10,height=10)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,bw=T,textSize=22,trans=F,ylab="% normalised reads")#+theme(text=element_text(size=28))+scale_y_continuous(labels=scaleFUN)
dev.off()


### core

## bacteria

mybiom <- mybiom_sb
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
t1 <- aggregate(t(otu_table(myfiltbiom )),by=list(sample_data(myfiltbiom )[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>0,myfiltbiom )
tax_table(myfiltbiom)[,15] <- sub("\\(g\\)","",tax_table(myfiltbiom)[,15])
tax_table(myfiltbiom) <- sub("_"," ",tax_table(myfiltbiom))
postscript("b_y2_core_final3.3.eps",width=10,height=10)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,bw=T,textSize=22,trans=T,margins=c(0.5,0.2,0.2,2) ,ylab="% variance stabilised reads")
dev.off()

## fungi

mybiom <- mybiom_sf
sample_data(mybiom)$condition <- sub("-.*","",sample_data(mybiom)$condition)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
t1 <- aggregate(t(otu_table(myfiltbiom )),by=list(sample_data(myfiltbiom )[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>0,myfiltbiom )
tax_table(myfiltbiom)[,15] <- sub("\\(g\\)","",tax_table(myfiltbiom)[,15])
tax_table(myfiltbiom) <- sub("_fam_I.*\\(","(",tax_table(myfiltbiom))
postscript("f_y2_core_final3.3.eps",width=10,height=10)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,bw=T,textSize=22,trans=T,margins=c(0.5,0.2,0.2,1.5) ,ylab="% variance stabilised reads")#+theme(text=element_text(size=28))+scale_y_continuous(labels=scaleFUN)
dev.off()



pdf("poster_bacteria_core.pdf",width=16,height=16)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,coloured=F,textSize=48,trans=T,margins=c(0.2,0.2,0.2,3))
dev.off()


pdf("poster_fungi_core.pdf",width=18,height=16)
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,coloured=F,textSize=48,trans=T)
dev.off()


#===============================================================================
#       taxonomy
#===============================================================================

### bacteria

## 2015


mybiom <- mybiom_fb
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

# percentage top 20 by abundance OTUs
otu_counts <- rowSums(otu_table(myfiltbiom))
otu_counts <- otu_counts[order(otu_counts,decreasing=T)]
sum(otu_counts[1:20])/(sum(otu_counts))

# OTU counts per phylum
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"phylum")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x[-1]

# Top 20 by genera
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ret_data=T,trans=F)


## 2016

mybiom <- mybiom_sb
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)$condition=="N",myfiltbiom)

# percentage top 20 by abundance OTUs
otu_counts <- rowSums(otu_table(myfiltbiom))
otu_counts <- otu_counts[order(otu_counts,decreasing=T)]
sum(otu_counts[1:20])/(sum(otu_counts))

# OTU counts per phylum
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"phylum")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x[-1]

# Top 20 by genera
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ret_data=T,trans=F)

# OTU counts per genus
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"genus")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x <-x[-1]
x$tots <- rowSums(x)
colnames(genus_sb)[1] <- "tots"




### fungi

## 2015


mybiom <- mybiom_ff
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

# percentage top 20 by abundance OTUs
otu_counts <- rowSums(otu_table(myfiltbiom))
otu_counts <- otu_counts[order(otu_counts,decreasing=T)]
sum(otu_counts[1:20])/(sum(otu_counts))

# top OTU counts
head(prop.table(otu_counts))


# OTU counts per phylum
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"phylum")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x[-1]

# Top 20 by genera
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ret_data=T,trans=F)


## 2016


mybiom <- mybiom_sf
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)$condition=="N",myfiltbiom)

# percentage top 20 by abundance OTUs
otu_counts <- rowSums(otu_table(myfiltbiom))
otu_counts <- otu_counts[order(otu_counts,decreasing=T)]
sum(otu_counts[1:20])/(sum(otu_counts))

head(prop.table(otu_counts))
    OTU453      OTU35       OTU2       OTU1     OTU594       OTU9
0.20966136 0.14776463 0.09551665 0.07273207 0.04581055 0.02433643


# OTU counts per phylum
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"phylum")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x[-1]

# Top 20 by genera
plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ret_data=T,trans=F)


#===============================================================================
#      Rank correlation
#===============================================================================

### bacteria

mybiom <- mybiom_fb
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

# OTU counts per genus
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"genus")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x <-x[-1]
x$tots <- rowSums(x)
genus_fb <- x[order(x$tots,decreasing=T),]

mybiom <- mybiom_sb
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

# OTU counts per genus
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"genus")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x <-x[-1]
x$tots <- rowSums(x)
genus_sb <- x[order(x$tots,decreasing=T),]


# rank correlation

genus_fb$genus <- rownames(genus_fb)
genus_sb$genus <- rownames(genus_sb)
test <- full_join(genus_fb[,5:6],genus_sb,by="genus")
test <- test[complete.cases(test),]

test <- test[,c(2,1,3)]

test <- test[order(test[,2],decreasing=T),]
test$fb_ord <- seq(1,nrow(test))

test <- test[order(test[,3],decreasing=T),]
test$sb_ord <- seq(1,nrow(test))

cor.test(~fb_ord+sb_ord,test,method="spearman")

### fungi

mybiom <- mybiom_ff
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

# OTU counts per genus
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"genus")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x <-x[-1]
x$tots <- rowSums(x)
genus_ff <- x[order(x$tots,decreasing=T),]

mybiom <- mybiom_sf
dds <- phylo_to_des(mybiom)
mybiom@otu_table@.Data <- counts(dds,normalize=T)
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

# OTU counts per genus
x <- sumTaxa(phylo_to_ubiom(myfiltbiom),"genus")
x[,2:ncol(x)] <- as.numeric(sapply(x[,2:ncol(x)],scaleFUN))
rownames(x) <- x[,1]
x <-x[-1]
x$tots <- rowSums(x)
genus_sf <- x[order(x$tots,decreasing=T),]


# rank correlation

genus_ff$genus <- rownames(genus_ff)
genus_sf$genus <- rownames(genus_sf)
test <- full_join(genus_ff[,5:6],genus_sf,by="genus")
test <- test[complete.cases(test),]

test <- test[,c(2,1,3)]

test <- test[order(test[,2],decreasing=T),]
test$ff_ord <- seq(1,nrow(test))

test <- test[order(test[,3],decreasing=T),]
test$sf_ord <- seq(1,nrow(test))

cor.test(~ff_ord+sf_ord,test,method="spearman")


#===============================================================================
#       Junk
#===============================================================================


myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
df <- plotTaxa(myfiltbiom,"rank","all",type=2, others=F,topn=20,ordered=T,coloured=F, ret_data=T)


mybiom <- mybiom_ff 
myfiltbiom <- prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)

t1 <- aggregate(t(otu_table(myfiltbiom )),by=list(sample_data(myfiltbiom )[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>0,myfiltbiom )

sample_data(myfiltbiom)$all <- "all"

tax_table(myfiltbiom) <- sub("_genera_"," ",tax_table(myfiltbiom))
tax_table(myfiltbiom) <- sub("_fam_I"," i",tax_table(myfiltbiom))
tax_table(myfiltbiom) <- sub("_"," ",tax_table(myfiltbiom))
tax_table(myfiltbiom) <- sub("_"," ",tax_table(myfiltbiom))


g <-  plotTaxa(myfiltbiom,"genus","all",type=2, others=F,topn=20,fitType="local",ordered=T,coloured=F)
g<-g+theme(text=element_text(size=18))
g <- g + theme(legend.position="none")
g
dev.off()


mycorebiom <- prune_taxa(apply(otu_prop_table,1,function(x) (sum(x>=min_freq))/ncol(otu_prop_table)>=min_samp),myfiltbiom)
