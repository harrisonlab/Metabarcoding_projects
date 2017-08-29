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
library(ggplot2)
library(reshape2)
library(tibble)

#===============================================================================
#       Load data 
#===============================================================================

# load denoised otu count table
countData <- read.table("final_16s_Rhiz_otu_table.txt",header=T,sep="\t",row.names=1, comment.char = "")

# load sample metadata
colData <- read.table("16S cultivar meta data.csv",header=T,sep=",",row.names=1)

# colnames are no the same in Rhiz and Endo data sets - altered here
rownames(colData) <- sub("E","R",rownames(colData))

# load taxonomy data
taxData <- read.table("16S_Rhiz.taxa",header=F,sep=",",row.names=1)

# reorder columns
taxData<-taxData[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]

# add best "rank" at 0.8 confidence and tidy-up the table
taxData<-phyloTaxaTidy(taxData,0.8)

# read phylip distance matrix 
mdat <- fread.phylip("16S_Rhiz.phy",skip=1)

# save data into a list, then ubiom_16S$countData to access countData and etc.
ubiom_BAC_Rhiz <- list(
	countData=countData,
	colData=colData,
	taxData=taxData,
	mdat=mdat,
	RHB="BAC_Rhiz"
)

# or all in one
ubiom_BAC_Endo <- list(
	countData=read.table("final_16s_Endo_otu_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("16S cultivar meta data.csv",header=T,sep=",",row.names=1),
	taxData=phyloTaxaTidy(read.table("16S_Rhiz.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.8),
	mdat=fread.phylip("16S_Endo.phy",skip=1),
	RHB="BAC_Endo"
) 

ubiom_FUN <- list(
	countData=read.table("final_ITs_otu_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
	colData=read.table("ITS cultivar meta data.csv",header=T,sep=",",row.names=1),
	taxData=phyloTaxaTidy(read.table("ITS.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.8),
	mdat=fread.phylip("ITS.phy",skip=1),
	RHB="FUN"
) 

# o.k. need to subset the fungal data  - therefore first attach all fungal objects 
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))

# construct new fungal endophyte list
ubiom_FUN_Endo <- list(
	countData=countData[,grep("E",colnames(countData))],
	colData=colData[grep("E",rownames(colData)),2, drop = FALSE],
	taxData=taxData,
	mdat=mdat,
	RHB="FUN_Endo"
)

# construct new fungal rhizoshere list
ubiom_FUN_Rhiz <- list(
	countData=countData[,grep("R",colnames(countData))],
	colData=colData[grep("R",rownames(colData)),2, drop = FALSE],
	taxData=taxData,
	mdat=mdat,
	RHB="FUN_Rhiz"
)

#===============================================================================
#       Create DEseq object
#===============================================================================

# attach objects - run only one
invisible(mapply(assign, names(ubiom_BAC_Rhiz), ubiom_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_BAC_Endo), ubiom_BAC_Endo, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_FUN_Endo), ubiom_FUN_Endo, MoreArgs=list(envir = globalenv())))
invisible(mapply(assign, names(ubiom_FUN_Rhiz), ubiom_FUN_Rhiz, MoreArgs=list(envir = globalenv())))

# ensure colData rows and countData columns have the same order
colData <- colData[colnames(countData),,drop = FALSE]

# simple Deseq design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# calculate size factors - use geoMeans function if every gene contains at least one zero 
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Read accumulation and filter
#===============================================================================

# calculate row sums from normalised read counts
df <- as.data.table(rowSums(counts(dds,normalized=T)),keep.rownames=F)

# calculate cumulative sum of ordered OTU counts 
#df <- cumsum(df[order(-V2)])
df <- cumsum(df[order(V1,decreasing=T)])

# add column name
colnames(df) <- RHB

# set cutoffs 
cutoffs = c(0.8,0.9,0.99,0.999)

# get cutoff poisiotns
mylines <- data.table(RHB=sapply(cutoffs,function(i) {nrow(df) - sum(df >= max(df,na.rm=T)*i)}),Cutoffs=cutoffs)

# log the count values
df <- log10(df)

# output pdf file
pdf(paste(RHB,"OTU_counts.pdf",sep="_"))

# create an empty ggplot object from the data table
g <- ggplot(data=df,aes_string(x=seq(1,nrow(df)),y=RHB))

# remove plot background and etc.
g <- g + theme_classic_thin()

# a colour pallete that can be descriminated by the colur blind 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot cumulative reads
g <- g + geom_line(size=1.5) + scale_colour_manual(values=cbbPalette) # single line so colour scale will only use first item

# add cutoff lines
g <- g + geom_vline(data=mylines,aes(xintercept=RHB),colour=cbbPalette[2:5])

# label cutoff lines
g <- g + geom_text(data = mylines, aes(label=Cutoffs,x=RHB, y = -Inf),angle = 90, inherit.aes = F, hjust = -.5, vjust = -.5)

# add axis lables
g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")

# print the plot
g

dev.off()

# o.k the above loses the OTU names - this keeps them, not worth editing the above.
df <- as.data.table(rowSums(counts(dds,normalize=T)),keep.rownames=T)
df <- df[order(-V2)]
df$V2 <- cumsum(df$V2)

# for Rhiz based on the graph we want the first 500 OTUs by abundance, Endo about 50 (or less)
myfilter <- df$V1[1:500] # bac rhiz
myfilter <- df$V1[1:50] # bac endo
myfilter <- df$V1[1:250] # fun endo
myfilter <- df$V1[1:150]  # fun rhiz

# filter out low abundance OTUs
dds <- dds[myfilter,]

# filter by proportion of reads
# myfilter <- (rowSums(counts(dds,normalize=T))/sum(counts(dds,normalize=T))) >= 0.00001

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
apply(mypca$x[,1:4],2,function(x) summary(aov(x~Cultivar,as.data.frame(dds@colData))))

# get sum of squares for all PC scores
sum_squares <- t(apply(mypca$x,2,function(x)t(summary(aov(x~Cultivar,as.data.frame(dds@colData)))[[1]][2])))

# name sum_squares columns
colnames(sum_squares) <- c("Cultivar","residual")

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
plotOrd(df,dds@colData,design="Cultivar",xlabel="PC1",ylabel="PC2")

# plot PC2 vs PC3
plotOrd(df,dds@colData,design="Cultivar",xlabel="PC2",ylabel="PC3",dimx="PC2",dimy="PC3")

# write to file
dev.off()

#===============================================================================
#       Diversity analysis
#===============================================================================

# ensure dds is filtered appropriately
dds <- dds[myfilter,]

#### ALPHA DIVERSITY ####

# function for converting normalised counts to integers and ceiling values below 1
alpha_counts <- function(X) {
	X[X==0]<-NA
	X[X<1] <- 1
	X[is.na(X)] <- 0
	return(round(X,0))
}

# transposed normalised read counts - converted to integers
OTU <- t(alpha_counts(counts(dds,normalize=T)))

# get alpha diversity measures with Vegan(Chao1, ACE, Shannon, Simpson)
all_alpha <- data.table(
	t(estimateR(OTU)),
	shannon = diversity(OTU, index = "shannon"),
	simpson = diversity(OTU, index = "simpson"),
	keep.rownames="Samples"
)

# rename samples to Cultivar name
all_alpha$Samples <- sub("[RE][0-9]$","",all_alpha$Samples)

# alpha diversity significance tests
sink(paste(RHB,"alpha.txt",sep="_"))
cat("summary(aov(S.chao1~Samples,all_alpha))","\n")
summary(aov(S.chao1~Samples,all_alpha))
cat("\n\nsummary(aov(S.ACE~Samples,all_alpha))","\n")
summary(aov(S.ACE~Samples,all_alpha))
cat("\n\nsummary(aov(shannon~Samples,all_alpha))","\n")
summary(aov(shannon~Samples,all_alpha))
cat("\n\nsummary(aov(simpson~Samples,all_alpha))","\n")
summary(aov(simpson~Samples,all_alpha))
sink()

### plot Alpha diversity

# get column names of all_alpha
measures = colnames(all_alpha)

# get standard error columns (Chao1 and ACE)
ses = colnames(all_alpha)[grep("^se\\.", colnames(all_alpha))]

# set measures to (measures - se) columns
measures = measures[!measures %in% ses]

# melt the data table by Sample
mdf <- melt(all_alpha,id.vars="Samples")

# add a column for se values
mdf$se <- as.double()

# temp holder for se labels
selabs = ses

# rename se labels to the same as the corresponding measures column
names(selabs) <- sub("se","S",selabs)

# add se values to se column for Chao1 and ACE
#mdf[variable %like% x]$se <- (sapply(names(selabs),function(x)mdf[variable %like% selabs[x]]$value))
invisible(sapply(names(selabs),function(x) mdf[variable %like% x]$se <<-  mdf[variable %like% selabs[x]]$value))

# remove se rows
mdf <- mdf[as.character(mdf$variable) %in% measures,]

# drop the levels (se column names) we've just removed 
mdf$variable <- droplevels(mdf$variable)

# create ggplot object
g <- ggplot(data=mdf,aes(x=Samples,y=value,colour=Samples))

# add a theme
g <- g + theme_classic_thin() %+replace% theme(	
	panel.border = element_rect(colour = "black", fill=NA, size=0.5),
	axis.text.x = element_text(angle = -90, vjust = 0.5,hjust = 0)
)

# add points
g <- g + geom_point(na.rm = TRUE,position = position_dodge(width = 0.5),size=1)

# add error bars
g <- g + geom_errorbar(aes(ymax = value + se, ymin = value -  se), width = 0.5,position = position_dodge(width = 0.5))

# add y label
g <- g + ylab("Alpha Diversity Measure")

# Add heading to each graph
g <- g + facet_wrap(~variable, nrow = 1,scales="free_y")

# output file
pdf(paste(RHB,"alpha.pdf",sep="_"))

# print the graphs
sapply(levels(mdf$variable),function(lev) print(g %+% mdf[variable==lev,]) )

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
adonis(dw~Cultivar,as.data.frame(dds@colData),parallel=12,permutations=9999,method="euclidean")
adonis(du~Cultivar,as.data.frame(dds@colData),parallel=12,permutations=9999,method="euclidean")
adonis(d5~Cultivar,as.data.frame(dds@colData),parallel=12,permutations=9999,method="euclidean")
sink()

### Plot beta diversity

# order rows and columns
d5 <- d5[sort(rownames(d5)),sort(colnames(d5))]

# output file
pdf(paste(RHB,"beta.pdf",sep="_"),height=7,width=7)

# heatmap plotting function
plotHeatmap(d5,textSize=11,axis.labels=T)

# close output file
dev.off()

#===============================================================================
#       Statistical analysis
#===============================================================================

# ensure dds is filtered appropriately
dds <- dds[myfilter,]

# p value for FDR cutoff
alpha <- 0.1

# drop any unused levels in Cultivar
dds$Cultivar <- droplevels(dds$Cultivar)

# add model to the DES object
design(dds) <- ~Cultivar

# calculate fit
dds <- DESeq(dds,parallel=T)

# prepare contrasts - down and dirty nested loop method, not very R
contrast <- ""
for (i in 1:length(levels(dds$Cultivar)[-1])){
	for( k in (i+1):length(levels(dds$Cultivar))) {
		x <- c("Cultivar",levels(dds$Cultivar)[i],levels(dds$Cultivar)[k])
		contrast<- rbind(contrast,x)
	}
}
contrast <- contrast[-1,]

# calculate results (this will take a while to run)
res.1.10 <- apply(contrast[1:10,],1,function(s) results(dds,alpha=alpha,parallel=T,contrast=s))

# add name to each result set (A vs B)
names(res.1.10) <- paste(contrast[1:10,2], contrast[1:10,3],sep=" vs ")

# extract 
logfc <- sapply(res,function(obj) data.frame(names(obj)=obj[,2]))
padj <- sapply(res,function(obj) data.frame(names(obj)=obj[,5]))
baseMean <- res[[1]][,1,drop = FALSE]

res.merge <- data.table(inner_join(data.table(OTU=rownames(baseMean),as.data.frame(baseMean)),data.table(OTU=rownames(taxData),taxData)))


lapply(res,function(obj) write.table(obj,paste(RHB,names(obj),"results.txt"
