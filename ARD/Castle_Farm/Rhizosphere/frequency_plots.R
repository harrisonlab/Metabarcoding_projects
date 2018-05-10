### Run after loading data from analysis_new.R ###

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

levels(md1$variable)[1] <- "Normalised"
levels(md2$variable)[1] <- "qPCR"

md <- rbind(md1,md2)
md$value <- as.numeric(md$value)
md$Samples <- md$variable
md$phylum <- sub(".*_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)

g <- ggplot(md[md$value>0.3,],aes(x=phylum,y=value,fill=Samples)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
#g <- g + geom_bar(stat="identity",colour="white",position = "dodge")
g <- g + scale_fill_manual(values = c("Normalised" = "black", "qPCR" = "orange"))
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

levels(md1$variable)[1] <- "Normalised"
levels(md2$variable)[1] <- "qPCR"

md <- rbind(md1,md2)
md$value <- as.numeric(md$value)
md$Samples <- md$variable
md$phylum <- sub("_"," ",md$phylum)
md$phylum  <- factor(md$phylum, levels=c(unique(md$phylum[md$phylum!="others"][order(md$value[md$phylum!="others"],decreasing=T)]),"others"),ordered=T)

g <- ggplot(md[md$value>0.3,],aes(x=phylum,y=value,fill=Samples)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
#g <- g + geom_bar(stat="identity",colour="white",position = "dodge")
g <- g + scale_fill_manual(values = c("Normalised" = "black", "qPCR" = "orange"))
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
#ggsave("FUN_frequency.pdf",g,width=8,height=7)
g2 <- g

### PUTTING THE PLOTS TOGETHER ###

g1 + theme(legend.position="hidden")
g2 + theme(legend.position="bottom")
library(gridExtra)
ggsave("Figure_2.pdf",grid.arrange(g1 + ggtitle("A")+ theme(legend.position="hidden") ,g2 + ggtitle("B")+ theme(legend.position="bottom"),nrow=2,padding = unit(-1, "line")),width=7,height=8)
