### Run after loading data from analysis_new.R ###

## FIGURE S2 ###
# number of taxa identified correctly at given confidence
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
dd <- sumTaxa(list(data.frame(cbind(taxData[,1,drop=F],1)[,2,drop=F]),taxData,data.frame(all=1)),conf=0.9,proportional=T,taxon="phylum")
colnames(dd)[2] <- "all"
dd <- rbind(dd[dd$all>=1,],c("others",sum(dd[dd$all<1,2])))
md <- melt(dd,id=colnames(dd)[1])
#levels(md$variable)[1] <- "OTUs"
md$value <- as.numeric(md$value)
#md$Normalisation <- md$variable
md$phylum <- sub(".*_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)
g <- ggplot(md,aes(x=phylum,y=value))
g <- g + geom_bar(stat="identity",colour="black")
g <- g  + xlab("")+ ylab("")
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
	plot.title = element_text(hjust = -0.11),
	axis.title.y=element_text(size=(14-2)))
ggsave("BAC_OTU_frequency.pdf",g,width=8,height=7)
g1 <- g

invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
dd <- sumTaxa(list(data.frame(cbind(taxData[,1,drop=F],1)[,2,drop=F]),taxData,data.frame(all=1)),conf=0.9,proportional=T,taxon="phylum")
colnames(dd)[2] <- "all"
dd <- rbind(dd[dd$all>=1,],c("others",sum(dd[dd$all<1,2])))
md <- melt(dd,id=colnames(dd)[1])
#levels(md$variable)[1] <- "OTUs"
md$value <- as.numeric(md$value)
#md$Normalisation <- md$variable
md$phylum <- sub(".*_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)
g <- ggplot(md,aes(x=phylum,y=value))
g <- g + geom_bar(stat="identity",colour="black")
g <- g  + xlab("")+ ylab("")
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
	plot.title = element_text(hjust = -0.11),
	axis.title.y=element_text(size=(14-2)))
ggsave("FUN_OTU_frequency.pdf",g,width=8,height=7)
g2 <- g

ggsave("Figure_S2.pdf",grid.arrange(g1,g2,left=textGrob(label="Frequency (%)",rot=90),nrow=2,ncol=1),width=7,height=9)


invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
sumTaxa(list(data.frame(cbind(taxData[,1,drop=F],1)[,2,drop=F]),taxData,data.frame(a=1)),conf=0.9,proportional=T,taxon="phylum")
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
taxData <- taxData[myfilter,]

sumTaxa(list(data.frame(cbind(taxData[,1,drop=F],1)[,2,drop=F]),taxData,data.frame(a=1)),conf=0.9,proportional=T,taxon="phylum")
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
sumTaxa(list(data.frame(cbind(taxData[,1,drop=F],1)[,2,drop=F]),taxData,data.frame(a=1)),conf=0.9,proportional=T,taxon="phylum")


library(gridExtra)
# bacteria
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
design<-~1

dds<-DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
nrm <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.9,proportional=T)
nrm <- rbind(nrm[nrm$all>=1,],c("others",sum(nrm[nrm$all<1,2])))

dds<-DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- 1/colData$bacq
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
qPCR <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.9,proportional=T)
qPCR <- rbind(qPCR[qPCR$all>=1,],c("others",sum(qPCR[qPCR$all<1,2])))

md1 <- melt(nrm,id=colnames(nrm)[1])
md2 <- melt(qPCR,id=colnames(qPCR)[1])

levels(md1$variable)[1] <- "DESeq2"
levels(md2$variable)[1] <- "qPCR"

md <- rbind(md1,md2)
md$value <- as.numeric(md$value)
md$Normalisation <- md$variable
md$phylum <- sub(".*_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)

g <- ggplot(md[md$value>0.3,],aes(x=phylum,y=value,fill=Normalisation)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
#g <- g + geom_bar(stat="identity",colour="white",position = "dodge")
g <- g + scale_fill_manual(values = c("DESeq2" = "black", "qPCR" = "orange","OTUs"="#56B4E9"))
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
	plot.title = element_text(hjust = -0.11),
	axis.title.y=element_text(size=(14-2)))
ggsave("BAC_frequency.pdf",g,width=8,height=7)
g1 <- g

# fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
design<-~1
dds<-DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
nrm <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.9,proportional=T)
nrm <- rbind(nrm[nrm$all>=1,],c("others",sum(nrm[nrm$all<1,2])))

dds<-DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- 1/colData$bacq
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
qPCR <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.9,proportional=T)
qPCR <- rbind(qPCR[qPCR$all>=1,],c("others",sum(qPCR[qPCR$all<1,2])))

md1 <- melt(nrm,id=colnames(nrm)[1])
md2 <- melt(qPCR,id=colnames(qPCR)[1])

levels(md1$variable)[1] <- "DESeq2"
levels(md2$variable)[1] <- "qPCR"

md <- rbind(md1,md2)
md$value <- as.numeric(md$value)
md$Normalisation <- md$variable
md$phylum <- sub("_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)

g <- ggplot(md[md$value>0.3,],aes(x=phylum,y=value,fill=Normalisation)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
#g <- g + geom_bar(stat="identity",colour="white",position = "dodge")
g <- g + scale_fill_manual(values = c("DESeq2" = "black", "qPCR" = "orange"))
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
	plot.title = element_text(hjust = -0.11),
	axis.title.y=element_text(size=(14-2)))
ggsave("FUN_frequency.pdf",g,width=8,height=7)
g2 <- g

### PUTTING THE PLOTS TOGETHER ###

ggsave("Figure_2.pdf",grid.arrange(g1 + ggtitle("A")+ theme(legend.position="hidden") ,g2 + ggtitle("B")+ theme(legend.position="bottom"),nrow=2,padding = unit(-1, "line")),width=7,height=8)



### PCA/NMDS figures ###
# bacteria - deseq normalised
mypca <- des_to_pca(dds)
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

g1 <- plotOrd(d,colData,design="Condition",axes=c(1,3),plot="label",label="Pair",labelSize=3.5,cbPalette=T,legend="hidden",ylims=c(-4,4))
g1 <- g1 + ggtitle("A")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(hjust = -0.08,size=14),legend.position="none",)


myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))

# add tree to phyloseq object
phy_tree(myphylo) <- njtree

# calculate NMDS ordination using weighted unifrac scores
set.seed(sum(utf8ToInt("Xiangming Xu")))
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)

g <- plotOrd(ordu$points,colData,design="Condition",plot="label",label="Pair",labelSize=4,cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",xlims=c(-0.1,0.1))
g$layers[[1]] <- NULL
g  <- g + geom_point(size = 0, stroke = 0)  # OR  geom_point(shape = "") +
g  <- g + geom_label(show.legend = FALSE,size=3.5)
g  <- g + guides(colour = guide_legend(override.aes = list(size = 5, shape = c(utf8ToInt("H"), utf8ToInt("S")))))
g  <- g + scale_colour_manual(name = "Condition", breaks = c("H","S"), labels = c("",""),values=c("#000000", "#E69F00"))
g2 <- g + ggtitle("B")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(hjust = -0.22,,size=14),legend.position="bottom")

ggsave("Figure_5.pdf",grid.arrange(g1,g2,nrow=2),width=7,height=8)




# Figure S1


invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
g1 <- plotCummulativeReads(counts(dds,normalize=T))
sizeFactors(dds) <- 1/colData$bacq
g2 <- plotCummulativeReads(counts(dds,normalize=T))

invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <- 1/colData$bacq
g3 <- plotCummulativeReads(counts(dds,normalize=T))

invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
g4 <- plotCummulativeReads(counts(dds,normalize=T))

invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
g5 <- plotCummulativeReads(counts(dds,normalize=T))

g1 <- g1 + ggtitle("Bacteria DESeq") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g2 <- g2 + ggtitle("Bacteria qPCR") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g3 <- g3 + ggtitle("Fungi qPCR") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g4 <- g4 + ggtitle("Oomycete") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g5 <- g5 + ggtitle("Nematode") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))

ggsave("Figure_S1.pdf",grid.arrange(g1,g2,g3,g4,g5,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=3,ncol=2),width=7,height=9)

# Table S2
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]

# remove low count and control samples
myfilter <- colData$Condition!="C"

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
