# plot_taxa.r 
# obj (phloseq) must is a phyloseq object which must include taxonomy and sample data
# taxon (str) is the taxonomic level of interest
# condition (str) describes how the samples should be grouped (must be column of sample data)
# proportional (bool) whether the graph should use proportional or absolute values
# cutoff (double) for proportional graphs. Taxons below this value will be pooled into "other"
# topn (int)taxons to display (by total reads) for non-prortional graphs. Taxons below topn will be pooled into "other"
# type is limited to by sample (1) or by taxa (2)
# fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
# ncol is a ggplot paramter to use n columns for the legend

plot_taxa <- function(obj,taxon="phylum",condition="condition",proportional=T,cutoff=1,topn=20,type=1,fixed=F,ncol=1) {
	suppressPackageStartupMessages(require(DESeq2))
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(scales))

	countData=as.data.frame(as.matrix(obj@otu_table)@.Data)
	taxa=obj@tax_table
	colData=as.data.frame(sample_data(obj)@.Data,row.names=sample_data(obj)@row.names)	
	colnames(colData) <- sample_data(obj)@names
	colnames(taxa) <- c("kingdom","phylum","class","order","family","genus","species")
	taxa <- sub("*._+","",taxa)

	dds <- 	DESeqDataSetFromMatrix(countData,colData,~1)
	if (sum(apply(countData,1,function(x) prod(x!=0)))>0) {
		suppressPackageStartupMessages(require(edgeR))
		sizeFactors(dds) <- calcNormFactors(counts(dds))	
	} else {
		print("every gene contains at least one zero")
		print("ignoring all zero values")
		gm_mean = function(x, na.rm=TRUE){
			exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
		}
		geoMeans = apply(counts(dds), 1, gm_mean)
		dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
	}
	countData <- as.data.frame(assay(varianceStabilizingTransformation(dds,blind=F,fitType="local")))
	countData[countData<0] <- 0

	mybiom <- list(
		countData=countData,
		taxa=taxa,
		colData=colData
	)


	taxa_sum <- sumTaxa(mybiom,taxon=taxon,condition=condition)
	###proportional
	if(proportional) {
		mybiom$colData$MLUflop <- 1 #assigns the MLU flop digit
		tx <- sumTaxa(mybiom,taxon=taxon,"MLUflop")
		tx[,-1] <- prop.table(as.matrix(tx[,-1]),2)*100
		txk <- tx[tx[,2]>=cutoff,1]
		taxa_sum[,-1] <- prop.table(as.matrix(taxa_sum[,-1]),2)*100
	} else {
		taxa_sum <- taxa_sum[order(rowSums(taxa_sum[,-1]),decreasing=T),]	
		txk <- taxa_sum[1:topn,1]
	}

	### common
	taxa_cut <- taxa_sum[taxa_sum[,1]%in%txk,]
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=T),]
	taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="other" ,t(colSums(taxa_sum[!taxa_sum[,1]%in%txk,-1]))),names(taxa_cut)))
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	md2 <- melt(taxa_cut,id=colnames(taxa_cut)[1])
	md2$variable <- factor(md2$variable, levels =levels(md2$variable)[order(levels(md2$variable))]  )

	if (type==1) {
		g <- ggplot(md2,aes_string(x=md2$variable,y=md2$value,fill=taxon))	
	} else if (type==2) {
		colnames(md2) <- c("taxa","all","value")
		g <- ggplot(md2,aes_string(x=as.factor(md2[,1]),y=md2[,3],fill="all"))
	} else {
	  g <- ggplot(md2,aes_string(x=md2$variable,y=md2$value,fill=taxon))
	}
	
	g <- g + geom_bar(stat="identity",colour="white")
	g <- g + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	g <- g  + xlab("")
	if (fixed) {
		g <- g  + coord_fixed(ratio = 0.1)
	} 
	g <- g + scale_y_continuous(expand = c(0, 0),labels = comma)
	g <- g + ylab("")
	g <- g + guides(fill=guide_legend(ncol=ncol))
	g <- g + theme_bw()
	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin=unit(c(0.2,0.2,0.2,1.5),"cm"), 
	    axis.line.y = element_line(colour = "black"),axis.ticks.x=element_blank())
	return(g)
}
