combine_biom <- function(locX,locY) {
	biom1 <- read.table(locX,header=T,sep="\t", comment.char="")	
	biom2 <- read.table(locY,header=T,sep="\t", comment.char="")
	biom <- merge(biom1,biom2,by.x="X.OTU.ID",by.y="X.OTU.ID",all=T)
	biom[is.na(biom)] <- 0
	return(biom)	
}

ddsCalc <- function(X, design=~condition) {
	suppressPackageStartupMessages(require(DESeq2))
	dds <- 	DESeqDataSetFromMatrix(X$countData,X$colData,design)
	if (sum(apply(X$countData,1,function(x) prod(x!=0)))>0) {
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
	return(DESeq(dds, fitType="local"))
} 

fltTaxon <- function(X,taxon="phylum") {
  n <- which(colnames(X$taxa)==taxon)
	x <- aggregate(X$countData,by=X$taxa[,1:n],sum)
	ls.biom <- list(x[,(n+1):ncol(x)],X$colData,x[,1:n])
	names(ls.biom) <- c("countData","colData","taxa")
	return(ls.biom)	
}

plotPCA <- function (	
			object, 
			intgroup = "condition",
			labelby,
			ntop = 500,
			pcx = 1,
			pcy = 2, 
			returnData = FALSE,
			cofix=F,
			transform= function(	object,
						blind=F,
						fitType="local"
					) 
			{
					suppressPackageStartupMessages(require(DESeq2))
					return(varianceStabilizingTransformation(object,blind=blind,fitType=fitType))
			}
		) 
{
    suppressPackageStartupMessages(require(genefilter))
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(DESeq2))
    rld <- transform(object=object)
    rv <- rowVars(assay(rld))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(rld)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(rld@colData))) {
        stop("the argument 'intgroup' should specify columns of colData")
    }
    intgroup.df <- as.data.frame(rld@colData[, intgroup,drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(rld)[[intgroup]]
    }

   shape <- if (missing(labelby)) {factor("All")}else {as.factor(rld@colData[,labelby]) }

    d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,intgroup.df,shape=shape)
    colnames(d)[grep("group", colnames(d))] <- intgroup
    colnames(d)[grep("shape", colnames(d))] <- labelby

    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }

    if(cofix) {
	d[,1] <- d[,1] * percentVar[pcx]
	d[,2] <- d[,2] * percentVar[pcy]
    }

    g <- ggplot()
    g <- g + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
    g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group, shape=shape),size=3)
    g <- g + scale_colour_discrete(name=intgroup)
    g <- g + scale_shape_discrete(name=labelby)
    g <- g + xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance"))
    g <- g + ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
    return(g)
}

plotPCAWithLabels <- function (object, intgroup = "condition", ntop = 500,pcx = 1,pcy = 2, returnData = FALSE)
{
    suppressPackageStartupMessages(require(genefilter))
    suppressPackageStartupMessages(require(ggplot2))
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(object@colData))) {
        stop("the argument 'intgroup' should specify columns of colData")
    }
    intgroup.df <- as.data.frame(object@colData[, intgroup,
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,
        intgroup.df, name = object$label)
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }

    ggplot() +
    geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=3) +
    geom_text(data=d, mapping=aes(x=PC1, y=PC2, label=name,colour=group), size=3, vjust=2, hjust=0.5) +
    xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance")) +
    ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
}

plotTaxa <- function(
			obj=mybiom, 	# obj (phloseq) must is a phyloseq object which must include taxonomy and sample data
			taxon="phylum", 	# taxon (str) is the taxonomic level of interest
			condition, 	# condition (str) describes how the samples should be grouped (must be column of sample data)
			proportional=T,	# proportional (bool) whether the graph should use proportional or absolute values
			cutoff=1, 	# cutoff (double) for proportional graphs. 
			topn=0, 		# topn (int)taxons to display (by total reads) for non-prortional graphs. T
			others=T, 	# combine values less than cutoff/topn into group "other"
			reorder=F, 	# order by value (max to min)
			type=1, 		# type is limited to by sample (1) or by taxa (2)
			fixed=F, 		# fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
			ncol=1, 		# ncol is a ggplot paramter to use n columns for the legend
			calcFactors="DES", # can accept a function e.g. to use edgeR size factors:
					   # function(counts){
					   #    suppressPackageStartupMessages(require(edgeR))
					   #    calcNormFactors(counts)
					   # }
			transform= function(	object,
						blind=F,
						fitType="local"
					) 
			{
					suppressPackageStartupMessages(require(DESeq2))
					return(varianceStabilizingTransformation(object,blind=blind,fitType=fitType))
			} 		# data transformation function 
		)
{
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
	if (class(calcFactors)!="function") {
		if (sum(apply(countData,1,function(x) prod(x!=0)))>0) {
			dds <- estimateSizeFactors(dds)
		} else {
			gm_mean = function(x, na.rm=TRUE){
				exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
			}
			geoMeans = apply(counts(dds), 1, gm_mean)
			dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
		}
	} else {
		sizeFactors(dds) <- calcFactors(counts(dds))
	}	

	countData <- as.data.frame(assay(transform(dds)))
	countData[countData<0] <- 0

	mybiom <- list(
		countData=countData,
		taxa=taxa,
		colData=colData
	)

	taxa_sum <- sumTaxa(mybiom,taxon=taxon,condition=condition)

	if(!topn) {
		mybiom$colData$MLUflop <- 1 #assigns the MLU flop digit
		tx <- sumTaxa(mybiom,taxon=taxon,"MLUflop")
		tx[,-1] <- prop.table(as.matrix(tx[,-1]),2)*100
		txk <- tx[tx[,2]>=cutoff,1]
	} else {
		taxa_sum[,ncol(taxa_sum)+1]<- 0
		taxa_sum <- taxa_sum[order(rowSums(taxa_sum[,-1]),decreasing=T),]
		taxa_sum <- taxa_sum[,-ncol(taxa_sum)]	
		txk <- taxa_sum[1:topn,1]
	}
	
	if(proportional) {
		taxa_sum[,-1] <- prop.table(as.matrix(taxa_sum[,-1]),2)*100
	}

	### common
	taxa_cut <- taxa_sum[taxa_sum[,1]%in%txk,]
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=T),]
	if(others) {
		taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="other" ,t(colSums(taxa_sum[!taxa_sum[,1]%in%txk,-1]))),names(taxa_cut)))
	}
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	if(reorder) {
		taxa_cut[,ncol(taxa_cut)+1] <- 0
		taxa_cut[,1] <- reorder(taxa_cut[,1],-rowSums(taxa_cut[,-1]))
		taxa_cut <- taxa_cut[,-ncol(taxa_cut)]
	}
	md2 <- melt(taxa_cut,id=colnames(taxa_cut)[1])
	md2$variable <- factor(md2$variable, levels =levels(md2$variable)[order(levels(md2$variable))]  )
	if (type==1) {
		g <- ggplot(md2,aes_string(x=md2[,2],y=md2[,3],fill=taxon))
	} else if (type==2) {
		colnames(md2) <- c("taxa","all","value")
		g <- ggplot(md2,aes_string(x=as.factor(md2[,1]),y=md2[,3],fill="all"))
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
	if (type==1) {
		g <- g + theme(legend.text=element_text(face="italic"))
	} else if (type==2) {
		g <- g + theme(axis.text.x=element_text(face="italic"))
	}	
	return(g)
}


phylo_to_des <- function(X, design=~condition) {
    suppressPackageStartupMessages(require(DESeq2))
    dds <-  phyloseq_to_deseq2(X,design)
    if (sum(apply(counts(dds),1,function(x) prod(x!=0)))>0) {
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
    return(DESeq(dds, fitType="local"))
} 

sumTaxa <- function(X,taxon="phylum",condition="condition") {
	suppressPackageStartupMessages(require(plyr))
	suppressPackageStartupMessages(require(reshape2))
	
	X$taxa <- as.data.frame(t(apply(X$taxa,1,function(x) {
			if (x[2]=="unknown") {x[2] <- paste(tx[1],"(k)",sep="")}
			if (x[3]=="unknown") {if(any(grep('\\(',x[2]))) {x[3]<-x[2]}else{x[3]<-paste(x[2],"(p)",sep="")}}
			if (x[4]=="unknown") {if(any(grep('\\(',x[3]))) {x[4]<-x[3]}else{x[4]<-paste(x[3],"(c)",sep="")}}
			if (x[5]=="unknown") {if(any(grep('\\(',x[4]))) {x[5]<-x[4]}else{x[5]<-paste(x[4],"(o)",sep="")}}
			if (x[6]=="unknown") {if(any(grep('\\(',x[5]))) {x[6]<-x[5]}else{x[6]<-paste(x[5],"(f)",sep="")}}
			if (x[7]=="unknown") {if(any(grep('\\(',x[6]))) {x[7]<-x[6]}else{x[7]<-paste(x[6],"(g)",sep="")}}
						
		return(x)}
	      )))
	tx <- X$taxa[,taxon]
	dtx <- cbind(X$countData,tx)
	md <- melt(dtx,id="tx")
	md$variable <- mapvalues(md$variable,from=rownames(X$colData), to=as.character(X$colData[,condition]))
	nd <- dcast(md,...~variable,sum)
	colnames(nd)[1] <- taxon
	return(nd)
}

ubiom <- function(locX,locY,locZ) {
	options(stringsAsFactors = FALSE)
	countData <- read.table(locX,header=T,sep="\t", comment.char="")
	rownames(countData ) <- countData [,1]
	countData <- countData [,-1]
	taxa <- read.csv(locY,header=F)
	taxa <- taxa [,c(1,2,4,6,8,10,12,14)]
	rownames(taxa) <- taxa[,1]
	taxa <- taxa[,-1]
	colnames(taxa) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
	colData <- read.table(locZ,sep="\t",header=T)
	rownames(colData) <- colData [,1]
	colData <- colData[,-1,drop=FALSE]
	countData <- countData[,rownames(colData)]
	ls.biom <- list(countData,colData, taxa)
	names(ls.biom) <- c("countData","colData","taxa")
	return(ls.biom)
}

vst <- function(object,blind=F,fitType="local") {
	suppressPackageStartupMessages(require(DESeq2))
	return(varianceStabilizingTransformation(object,blind=blind,fitType=fitType))
}


