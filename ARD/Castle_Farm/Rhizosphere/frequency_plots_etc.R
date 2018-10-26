### Run after loading data from analysis_new.R ###

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
sizeFactors(dds) <- colData$bacq
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
qPCR <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.9,proportional=T)
qPCR <- rbind(qPCR[qPCR$all>=1,],c("others",sum(qPCR[qPCR$all<1,2])))

md1 <- melt(nrm,id=colnames(nrm)[1])
md2 <- melt(qPCR,id=colnames(qPCR)[1])

levels(md1$variable)[1] <- "Raw"
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
sizeFactors(dds) <- colData$funq
dds <- dds[rowSums(counts(dds, normalize=T))>4,]
qPCR <- sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),conf=0.9,proportional=T)
qPCR <- rbind(qPCR[qPCR$all>=1,],c("others",sum(qPCR[qPCR$all<1,2])))

md1 <- melt(nrm,id=colnames(nrm)[1])
md2 <- melt(qPCR,id=colnames(qPCR)[1])

levels(md1$variable)[1] <- "Raw"
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



### PCA/NMDS single figures ###

# bacteria - deseq normalised
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
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))
phy_tree(myphylo) <- njtree
set.seed(sum(utf8ToInt("Xiangming Xu")))
ordu_bac_norm = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
mypca_bac_norm <- des_to_pca(dds)


mypca <- des_to_pca(dds)
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

g <- plotOrd(ordu_bac_norm$points,colData,design="Condition",plot="label",label="Pair",labelSize=4,cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",xlims=c(-0.08,0.08))
g1 <- g + ggtitle("A")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(hjust = -0.17,size=14),legend.position="none")

g <- plotOrd(d,colData,design="Condition",axes=c(1,3),plot="label",label="Pair",labelSize=3.5,cbPalette=T,legend="hidden",ylims=c(-4,4))
g$layers[[1]] <- NULL
g  <- g + geom_point(size = 0, stroke = 0)  # OR  geom_point(shape = "") +
g  <- g + geom_label(show.legend = FALSE,size=3.5)
g  <- g + guides(colour = guide_legend(override.aes = list(size = 5, shape = c(utf8ToInt("H"), utf8ToInt("S")))))
g  <- g + scale_colour_manual(name = "Condition", breaks = c("H","S"), labels = c("",""),values=c("#000000", "#E69F00"))
g2 <- g + ggtitle("B")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(hjust = -0.12,size=14),legend.position="bottom",)

ggsave("Figure_5BA.pdf",grid.arrange(g1,g2,nrow=2),width=7,height=8)



# Figure S2


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
sizeFactors(dds) <- colData$bacq
g2 <- plotCummulativeReads(counts(dds,normalize=T))

invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
g3 <- plotCummulativeReads(counts(dds,normalize=T))
sizeFactors(dds) <- colData$funq
g4 <- plotCummulativeReads(counts(dds,normalize=T))

invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]
g5 <- plotCummulativeReads(counts(dds,normalize=T))
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <- left_join(colData,ubiom_FUN$colData)$funq
myfilter <- row.names(countData[row.names(countData) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]
g6 <- plotCummulativeReads(counts(dds,normalize=T))
			    
			    
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
colData <- colData[names(countData),]
myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
exclude<-which(!myfilter)
myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
colData <- droplevels(colData[myfilter,])
countData <- countData[,myfilter]
dds<-DESeqDataSetFromMatrix(countData,colData,~1)
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
g7 <- plotCummulativeReads(counts(dds,normalize=T))

g1 <- g1 + ggtitle("Bacteria Raw") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g2 <- g2 + ggtitle("Bacteria qPCR") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g3 <- g3 + ggtitle("Fungi Raw") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g4 <- g4 + ggtitle("Fungi qPCR") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g5 <- g5 + ggtitle("Oomycete Raw") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g6 <- g6 + ggtitle("Oomycete qPCR") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))
g7 <- g7 + ggtitle("Nematode") + theme_classic_thin() %+replace% theme(axis.title=element_blank())#,plot.title = element_text(size=14))

ggsave("NEW_Figure_S2.pdf",grid.arrange(g1,g2,g3,g4,g5,g6,g7,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=4,ncol=2),width=8,height=10)

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


### NMDS/PCA combined plots ### 		
			    
qf <- function(colData,countData,sf="norm",rd=F) {
	colData <- colData[names(countData),]
	myfilter <- (colSums(countData)>=1000) & colData$Condition!="C"
	exclude<-which(!myfilter)
	myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)
	colData <- droplevels(colData[myfilter,])
	countData <- countData[,myfilter]
	design<-~1
	dds<-DESeqDataSetFromMatrix(countData,colData,design)
	if(sf=="norm"){
		sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))
	} else	{
		sizeFactors(dds) <- colData[[sf]]
	}
	dds <- dds[rowSums(counts(dds, normalize=T))>4,]
	if(rd) return(dds)
	myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))
	phy_tree(myphylo) <- njtree
	set.seed(sum(utf8ToInt("Xiangming Xu")))
	ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
	mypca <- des_to_pca(dds)
	return(list(mypca=mypca,ordu=ordu,colData=colData))
}


invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
bac_norm <- qf(colData,countData)
bac_qpcr <- qf(colData,countData,"bacq")

invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
fun_norm <- qf(colData,countData)
fun_qpcr <- qf(colData,countData,"funq")

invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
dds <- qf(colData,countData,rd=T)
myfilter <- row.names(counts(dds)[row.names(counts(dds)) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))
phy_tree(myphylo) <- njtree
set.seed(sum(utf8ToInt("Xiangming")))
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
mypca <- des_to_pca(dds)
oo_norm <- list(mypca=mypca,ordu=ordu,colData=colData(dds))

invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
dds <- qf(colData,countData,rd=T)
sizeFactors(dds) <- left_join(as.data.frame(colData(dds)),ubiom_FUN$colData)$funq
myfilter <- row.names(counts(dds)[row.names(counts(dds)) %in% row.names(taxData[(taxData$kingdom=="SAR"|as.numeric(taxData$k_conf)<=0.5),]),])
dds <- dds[myfilter,]
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))
phy_tree(myphylo) <- njtree
set.seed(sum(utf8ToInt("Xiangming Xu")))
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
mypca <- des_to_pca(dds)
oo_ooq <- list(mypca=mypca,ordu=ordu,colData=colData(dds))

invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
dds<- qf(colData,countData,rd=T)
myfilter <- row.names(taxData[as.number(taxData$c_conf)>0.9 & as.number(taxData$o_conf)>0.9,])
dds <- dds[rownames(dds)%in%myfilter,]
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),taxData,colData))
phy_tree(myphylo) <- njtree
set.seed(sum(utf8ToInt("Xiangming Xu")))
ordu = ordinate(myphylo, "NMDS", "unifrac", weighted=TRUE)
mypca <- des_to_pca(dds)
nem_norm <- list(mypca=mypca,ordu=ordu,colData=colData(dds))


#NMDS 
g <- plotOrd(bac_norm$ordu$points,bac_norm$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2")+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g1 <- g + ggtitle("Bacteria_raw")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

g <- plotOrd(bac_qpcr$ordu$points,bac_qpcr$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2")+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g2 <- g + ggtitle("Bacteria_qPCR")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

g <- plotOrd(fun_norm$ordu$points,fun_norm$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",ylims=c(-0.15,0.15))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g3 <- g + ggtitle("Fungi_raw")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

g <- plotOrd(fun_qpcr$ordu$points,fun_qpcr$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",xlims=c(-0.3,0.3),ylims=c(-0.15,0.15))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g4 <- g + ggtitle("Fungi_qPCR")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

g <- plotOrd(oo_norm$ordu$points,oo_norm$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",ylims=c(-0.08,0.08 ))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g5 <- g + ggtitle("Oomycete_raw")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

g <- plotOrd(oo_qpcr$ordu$points,oo_qpcr$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",ylims=c(-0.08,0.08 ))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g6 <- g + ggtitle("Oomycete_qPCR")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

g <- plotOrd(nem_norm$ordu$points,nem_norm$colData,design="Condition",facet="Pair",cbPalette=T,legend="bottom",xlabel="NMDS1",ylabel="NMDS2",ylims=c(-0.25,0.25))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g7 <- g + ggtitle("Nematode")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="right")
ggsave("Figure_S10_NEW.pdf",grid.arrange(g1,g2,g3,g4,g5,g6,g7,nrow=4),width=7,height=9)

# PCA
d <-t(data.frame(t(bac_norm$mypca$x)*bac_norm$mypca$percentVar))
g1 <- plotOrd(d,bac_norm$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden",ylims=c(-3.5,4.5))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g1 <- g1 + ggtitle("Bacteria_raw")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

d <-t(data.frame(t(bac_qpcr$mypca$x)*bac_qpcr$mypca$percentVar))
g2 <- plotOrd(d,bac_qpcr$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden",xlims=c(-15,15))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g2 <- g2 + ggtitle("Bacteria_qPCR")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

d <-t(data.frame(t(fun_norm$mypca$x)*fun_norm$mypca$percentVar))
g3 <- plotOrd(d,fun_norm$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden",ylims=c(-4.5,3.5))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g3 <- g3 + ggtitle("Fungi_raw")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

d <-t(data.frame(t(fun_qpcr$mypca$x)*fun_qpcr$mypca$percentVar))
g4 <- plotOrd(d,fun_qpcr$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden",ylims=c(-4,4))+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g4 <- g4 + ggtitle("Fungi_qPCR")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

d <-t(data.frame(t(oo_norm$mypca$x)*oo_norm$mypca$percentVar))
g5 <- plotOrd(d,oo_norm$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden")+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g5 <- g5 + ggtitle("Oomycete_raw")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

d <-t(data.frame(t(oo_qpcr$mypca$x)*oo_qpcr$mypca$percentVar))
g6 <- plotOrd(d,oo_qpcr$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden")+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g6 <- g6 + ggtitle("Oomycete_qPCR")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="none")

d <-t(data.frame(t(nem_norm$mypca$x)*nem_norm$mypca$percentVar))
g7 <- plotOrd(d,nem_norm$colData,design="Condition",axes=c(1,2),facet="Pair",cbPalette=T,legend="hidden")+ geom_line(aes(group=facet),alpha=0.125,linetype=3,colour="#000000")
g7 <- g7 + ggtitle("Nematode")+ theme_classic_thin(base_size=12)%+replace% theme(plot.title = element_text(size=14),legend.position="right")

ggsave("Figure_S11_NEW.pdf",grid.arrange(g1,g2,g3,g4,g5,g6,g7,nrow=4),width=9,height=10)




#===============================================================================
#       Sample rarefaction plots
#===============================================================================        

library(grid)
library(gridExtra)
library(viridis)
				    
gfunc <- function(countData,coldata,title) {        
  colData <- colData[names(countData),]

  # remove low count and control samples
  myfilter <- colData$Condition!="C"

  # remove Pair of any sample with a low count
  exclude<-which(!myfilter)
  myfilter <- myfilter&sapply(colData$Pair,function(x) length(which(x==colData$Pair[-exclude]))>1)

  # apply filter
  colData <- droplevels(colData[myfilter,])
  countData <- countData[,myfilter]

  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))

  # get cummulative sum of each sample
  DT <- cumsum(DT)    

  # log the count values                            
  DT <- log10(DT)

  # relabel columns
  colnames(DT) <- sub("(X)([0-9]+[HS])(.*)","\\2",colnames(DT))

  # set values larger than maximum for each column to NA
  DT <- data.table(apply(DT,2,function(x) {x[(which.max(x)+1):length(x)]<- NA;x}))
  
  # remove rows with all NA
  DT <- DT[rowSums(is.na(DT)) != ncol(DT), ]
  
  # add a count column to the data table
  DT$x <- seq(1,nrow(DT))
                             
  # melt the data table for easy plotting 
  MDT <- melt(DT,id.vars="x")
			      
  # create an empty ggplot object from the data table
  g <- ggplot(data=MDT,aes(x=x,y=value,colour=variable))

  # remove plot background and etc.
  g <- g + theme_classic_thin() %+replace% theme(legend.position="none",axis.title=element_blank())

  # plot cumulative reads
  g <- g + geom_line(size=1.5) + scale_colour_viridis(discrete=T)

  # add axis lables
  g <- g + ggtitle(title)
  #g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")

  # print the plot
  g
}                        

invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
g1 <- gfunc(countData,colData,"Bacteria")
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
g2 <- gfunc(countData,colData,"fungi")
invisible(mapply(assign, names(ubiom_OO), ubiom_OO, MoreArgs=list(envir = globalenv())))
g3 <- gfunc(countData,colData,"Oomycetes")
invisible(mapply(assign, names(ubiom_NEM), ubiom_NEM, MoreArgs=list(envir = globalenv())))
g4 <- gfunc(countData,colData,"Nematodes")

glegend <- get_legend(g)
ggsave("rarefaction_all.pdf",grid.arrange(g1,g2,g3,g4,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=2))                              

#===============================================================================
#       qPCR ratio plot
#===============================================================================   
# se method borrowed from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
test <- melt(colData[1:30,c(1,2,5)],idvars="Condition")
test <- test[!test$Pair==3,]
#colNames(test)[4] <- "Ratio"
mydf <- ddply(test, "Condition", .drop=.drop,.fun = function(x, col) {
  c(N    = length2(x[[col]], na.rm=na.rm),
    Ratio = mean   (x[[col]], na.rm=na.rm),
    sd   = sd     (x[[col]], na.rm=na.rm)
   )
}, measurevar)
mydf$se <- datac$sd / sqrt(datac$N) 

g <- ggplot(mydt,aes(x=Condition,y=Ratio,fill=Condition))
g <- g + geom_bar(position=position_dodge(), stat="identity") +
g <- g + geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.2, position=position_dodge(.9))
g <- g + theme_blank(base_size=12) %+replace%
g <- g + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5,linetype=1),
  axis.title.x=element_blank(),
  legend.position="none",
  axis.text.x = element_text(angle = 0, vjust = 1,hjust = 1),
  axis.ticks=element_blank())
# remove space at bottom of grapha and retain space at top
g <- g + scale_y_continuous(expand = expand_scale(mult = c(0, .075)))
ggsave("ratio.pdf",g)



   
