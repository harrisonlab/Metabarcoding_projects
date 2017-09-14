#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)
library(dplyr)
library(plyr)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions")
library(vegan)
library(ape)
library(GUniFrac)
library(gtools)
library(ggplot2)
library(tibble)

#===============================================================================
#       Load data 
#===============================================================================

# load denoised otu count table
countData <- read.table("final_16S_otu_table.txt",header=T,sep="\t",row.names=1, comment.char = "")

# load sample metadata
colData <- read.table("16S rotation meta data.csv",header=T,sep=",",row.names=1,colClasses=c("character","factor","factor"))

# load taxonomy data
taxData <- read.table("16S.taxa",header=F,sep=",",row.names=1)

# reorder columns
taxData<-taxData[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]

# add best "rank" at 0.8 confidence and tidy-up the table
taxData<-phyloTaxaTidy(taxData,0.8)

# read phylip distance matrix 
mdat <- fread.phylip("16S.phy",skip=1)

# save data into a list, then ubiom_16S$countData to access countData and etc.
ubiom_BAC <- list(
	countData=countData,
	colData=colData,
	taxData=taxData,
	mdat=mdat,
	RHB="BAC"
)

# add dds object to ubiom list
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC)

ubiom_FUN <- list(
	countData=read.table("final_ITs_otu_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("ITS rotation meta data.csv",header=T,sep=",",row.names=1,colClasses=c("character","factor","factor")),
	taxData=phyloTaxaTidy(read.table("ITS.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.8),
	mdat=fread.phylip("ITS.phy",skip=1),
	RHB="FUN"
)
ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN)

#===============================================================================
#       Attach objects
#===============================================================================

# attach objects - run only one
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Read accumulation and filter
#===============================================================================

# output pdf file
pdf(paste(RHB,"OTU_counts.pdf",sep="_"))

# plot cummulative reads
plotCummulativeReads(counts(dds,normalize=T))

# close pdf
dev.off()

#### Select filter ####
df <- as.data.table(rowSums(counts(dds,normalize=T)),keep.rownames=T)
df <- df[order(-V2)]

# from the cummulative plot select cut-off point
myfilter <- df$V1[1:1500] # bac
myfilter <- df$V1[1:500] # fun

# filter out low abundance OTUs
dds <- dds[myfilter,]

#===============================================================================
#       PCA
#===============================================================================

# perform PC decompossion on DES object
mypca <- des_to_pca(dds)

### ANOVA ###

# write to file
sink(paste(RHB,"PCA.txt"))

cat("
# ANOVA of first four PC scores \n")
apply(mypca$x[,1:4],2,function(x) summary(aov(x~Site + Year + Site * Year,as.data.frame(dds@colData))))

# get sum of squares for all PC scores
sum_squares <- t(apply(mypca$x,2,function(x)t(summary(aov(x~Site + Year + Site * Year,as.data.frame(dds@colData)))[[1]][2])))

# name sum_squares columns
colnames(sum_squares) <- c("Site","Year","Site:Year","residual")

# proportion of total sum of squares for PC
perVar <- t(apply(sum_squares,1,prop.table)) * mypca$percentVar

# check - should equal 1
# sum(colSums(perVar)) 

cat("
# % total variance explained by the aov model and residual \n")
colSums(perVar)/sum(colSums(perVar))*100

# end write to file
sink()

### Plot ###

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
df  <- t(data.frame(t(mypca$x)*mypca$percentVar))

# output pdf
pdf(paste(RHB,"PCA.pdf",sep="_"))

# plot PC1 vs PC2
plotOrd(df,dds@colData,design="Year",shape="Site",xlabel="PC1",ylabel="PC2",labels=T,textSize=8)

# plot PC2 vs PC3
plotOrd(df,dds@colData,design="Year",shape="Site",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3")

# write to file
dev.off()

# exclude variability "explained" site
pc.res <- resid(aov(mypca$x~Site,dds@colData))

# as above for residual values
d <- t(data.frame(t(pc.res)*mypca$percentVar))
		       
# plot PC1 vs PC2
plotOrd(d,dds@colData,design="Year",shape="Site",xlabel="PC1",ylabel="PC2",labels=F,textSize=8)
		       
# Bacteria has huge differences in year for PC1 for Sh site - obscures the same for the Ak site - worth splitting the data into two

# Fungal data has two big outliers Sh15ye and Sh20yc
# possibly due to large proportion of reads in a single OTU 
# X<-apply(counts(dds,normalize=T),2,prop.table)
# apply(X,2,max)[order(apply(X,2,max),decreasing=T)]

# for AK site
mypca <- des_to_pca(dds[,dds$Site=="Ak"])
df  <- t(data.frame(t(mypca$x)*mypca$percentVar))
pdf(paste(RHB,"Ak_PCA.pdf",sep="_"))
plotOrd(df,dds@colData[dds$Site=="Ak",],design="Year",xlabel="PC1",ylabel="PC2")
plotOrd(df,dds@colData[dds$Site=="Ak",],design="Year",xlabel="PC1",ylabel="PC2",cluster=0.80,centers=1)
plotOrd(df,dds@colData[dds$Site=="Ak",],design="Year",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3")
plotOrd(df,dds@colData[dds$Site=="Ak",],design="Year",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3",cluster=0.80,centers=1)
dev.off()

# for Sh site
mypca <- des_to_pca(dds[,dds$Site=="Sh"])
df  <- t(data.frame(t(mypca$x)*mypca$percentVar))

pdf(paste(RHB,"Sh_PCA.pdf",sep="_"))
plotOrd(df,dds@colData[dds$Site=="Sh",],design="Year",xlabel="PC1",ylabel="PC2")
plotOrd(df,dds@colData[dds$Site=="Sh",],design="Year",xlabel="PC1",ylabel="PC2",cluster=0.80,centers=1)
plotOrd(df,dds@colData[dds$Site=="Sh",],design="Year",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3")
plotOrd(df,dds@colData[dds$Site=="Sh",],design="Year",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3",cluster=0.80,centers=1)

dev.off()

#===============================================================================
#       Diversity analysis
#===============================================================================

# ensure dds is filtered appropriately
dds <- dds[myfilter,]

#### ALPHA DIVERSITY ####

all_alpha <- plot_alpha(counts(dds,normalize=T),returnData=T)

# join alpha diversity measures with sample meta data
all_alpha <- inner_join(all_alpha,as.data.table(as.data.frame(dds@colData),keep.rownames="Samples"))

# alpha diversity significance tests
sink(paste(RHB,"alpha.txt",sep="_"))
cat("summary(aov(S.chao1~Site+Year+Site*Year,all_alpha))","\n")
summary(aov(S.chao1~Samples,all_alpha))
cat("\n\nsummary(aov(S.ACE~Site+Year+Site*Year,all_alpha))","\n")
summary(aov(S.ACE~Samples,all_alpha))
cat("\n\nsummary(aov(shannon~Site+Year+Site*Year,all_alpha))","\n")
summary(aov(shannon~Samples,all_alpha))
cat("\n\nsummary(aov(simpson~Site+Year+Site*Year,all_alpha))","\n")
summary(aov(simpson~Site+Year+Site*Year,all_alpha))
sink()

### plot Alpha diversity

# output file
pdf(paste(RHB,"alpha.pdf",sep="_"))

# produce graphs
g <- plot_alpha(counts(dds,normalize=T),dds@colData,design="Site",colour="Year")

# print the graphs
sapply(levels(g$data$variable),function(lev) print(g %+% subset(g$data,variable==lev)))

dev.off()

#### BETA DIVERSITY (of unifrac distance) ####

# filter mdat 
mdat <- mdat[myfilter,myfilter]

# calculate rooted neighbour joining tree
mytree <- root(nj(as.dist(mdat)),outgroup=1,resolve=T)

# calculate unifrac distances for each sample
unifracs <- GUniFrac(t(counts(dds,normalize=T)),mytree,alpha=c(0, 0.5, 1))$unifracs

# extract weighted UniFrac distance
dw <- unifracs[, , "d_1"]

# extract unweighted UniFrac distance
du <- unifracs[, , "d_UW"] 

# Extract GUniFrac with alpha 0.5 "best method"
d5 <- unifracs[, , "d_0.5"]             

# test UniFracs for significane
sink(paste(RHB,"unifrac_adonis.txt",sep="_"))
adonis(dw~Site+Year+Site*Year,as.data.frame(dds@colData),parallel=12,permutations=9999,method="euclidean")
adonis(du~Site+Year+Site*Year,as.data.frame(dds@colData),parallel=12,permutations=9999,method="euclidean")
adonis(d5~Site+Year+Site*Year,as.data.frame(dds@colData),parallel=12,permutations=9999,method="euclidean")
sink()

### Plot beta diversity

# order rows and columns
#as.numeric(gsub("[A-Za-z]","",colnames(d5)))
d5 <- d5[mixedsort(rownames(d5)),mixedsort(colnames(d5))]

# output file
pdf(paste(RHB,"beta.pdf",sep="_"),height=7,width=7)

# heatmap plotting function - best to split into the two sites
plotHeatmap(d5[sub("[0-9].*","",rownames(d5))=="Ak",sub("[0-9].*","",colnames(d5))=="Ak"],textSize=11,axis.labels=T)
plotHeatmap(d5[sub("[0-9].*","",rownames(d5))=="Sh",sub("[0-9].*","",colnames(d5))=="Sh"],textSize=11,axis.labels=T)

# close output file
dev.off()

#===============================================================================
#       Statistical analysis
#===============================================================================


# ensure dds is filtered appropriately
dds <- dds[myfilter,]

# ensure Year is a factor
if(class(dds$Year)!="factor"){dds$Year <- as.factor$dds$Year}

# p value for FDR cutoff
alpha <- 0.1

# drop any unused levels in Cultivar
#dds$Cultivar <- droplevels(dds$)

# add model to the DES object - we're not interested at an individual site effect so don't add it to the model
design(dds) <- ~Site + Year #+ Site:Year

# calculate fit
dds <- DESeq(dds,parallel=T)

# prepare contrasts - down and dirty nested loop method, not very R
contrast <- ""
for (i in 1:length(levels(dds$Year)[-1])){
	for( k in (i+1):length(levels(dds$Year))) {
		x <- c("Year",levels(dds$Year)[i],levels(dds$Year)[k])
		contrast<- rbind(contrast,x)
	}
}
contrast <- contrast[-1,]

# calculate results (this will take a while to run)
res <- apply(contrast,1,function(s) results(dds,alpha=alpha,parallel=T,contrast=s))

# add name to each result set (A vs B)
names(res) <- paste(contrast[,2], contrast[,3],sep="_vs_")

# change names of log FC and padj column sto reflect the contrast
res <- Map(function(x, i) {
	colnames(x)[2]<-paste("Year",i, colnames(x)[2],sep="_");
	colnames(x)[6]<-paste("Year",i, colnames(x)[6],sep="_");
	return(x)
},res, names(res))

# merge contrast into a single table
m <- Reduce(function(...) inner_join(..., all=T),  lapply(res,function(obj) as.data.table(as.data.frame(obj[,c(1,2,6)]),keep.rownames="OTU")))

# merge with taxonomy data
res.merge <- data.table(inner_join(m,data.table(OTU=rownames(taxData),taxData)))

# write output
write.table(res.merge,paste(RHB,"diff_otu.txt",sep="_"),row.names=F,quote=F,na="",sep="\t")



