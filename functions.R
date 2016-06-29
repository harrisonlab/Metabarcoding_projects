
s3biom <- function(locX, locY) {
	require(biom)
	print("loading biom data")
	bm.biom <- read_biom(locX)
	colData <- read.table(locY)
	print("extracting useful data")
	df.biom.data <- data.frame(as.matrix(biom_data(bm.biom)))
	df.biom.data <- df.biom.data[,rownames(colData)]
	if (length(bm.biom$rows) > 7) {
		df.biom.taxa <-  suppressWarnings(as.data.frame(do.call(rbind,observation_metadata(bm.biom, 1:length(bm.biom$rows)))))
	} else {
		df.biom.taxa <-  suppressWarnings(as.data.frame(do.call(cbind,observation_metadata(bm.biom, 1:length(bm.biom$rows)))))
	}
	colnames(df.biom.taxa) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
	df.biom.data[is.na(df.biom.data)] <- c("No blast hit")
	ls.biom <- list(df.biom.data,colData, df.biom.taxa)
	names(ls.biom) <- c("countData","colData","taxa")
	return(ls.biom)
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

combine_biom <- function(locX,locY) {
	biom1 <- read.table(locX,header=T,sep="\t", comment.char="")	
	biom2 <- read.table(locY,header=T,sep="\t", comment.char="")
	biom <- merge(biom1,biom2,by.x="X.OTU.ID",by.y="X.OTU.ID",all=T)
	biom[is.na(biom)] <- 0
	return(biom)	
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

plotPCA <- function (object, intgroup = "condition", ntop = 500,pcx = 1,pcy = 2, returnData = FALSE)
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
    xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance")) +
    ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
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


fltTaxon <- function(X,taxon="phylum") {
  n <- which(colnames(X$taxa)==taxon)
	x <- aggregate(X$countData,by=X$taxa[,1:n],sum)
	ls.biom <- list(x[,(n+1):ncol(x)],X$colData,x[,1:n])
	names(ls.biom) <- c("countData","colData","taxa")
	return(ls.biom)	
}


##### TEST FUNCTIONS

#{
#	if (plotme) {
#		rld <- varianceStabilizingTransformation(dds,blind=F,fitType="local")
#		print("plotting PCA")
#		pdf(out,height=8,width=8)
#		plotPCAWithLabels(rld,)
#		dev.off()
#	}
#	)	
#}


ddsfun <- function(X,design=~condition) {
	suppressPackageStartupMessages(require(DESeq2))
	DESeqDataSetFromMatrix(X$countData,X$colData,design)
}

testfun <- function(X,plotme=F,out="funtest.pdf") {
#	if (plotme) {
		rld <- DESeqTransform(SummarizedExperiment(log2(counts(estimateSizeFactors(dds,geoMeans = geoMeans), normalized=TRUE) + 1),colData=colData(dds)))
		#varianceStabilizingTransformation(dds,blind=F,fitType="local")
		#pdf(out,height=8,width=8)
		plotPCAWithLabels(rld)
		#dev.off()
#	}
}

ddsMedian <- function(X, design=~condition, cutoff=3,plot=F,out="median.pdf") {
	suppressPackageStartupMessages(require(DESeq2))
	countData<- X$countData[rowSums(X$countData)>=cutoff,]
	dds <- 	DESeqDataSetFromMatrix(countData,X$colData,design)
	sizeFactors(dds) <- calcNormFactors(counts(dds))
	if (plot) {
		desPCAplotter(dds,out)
	}
	return(DESeq(dds, fitType="local"))
}


ddsGeoMeans <- function(X, design=~condition, cutoff=3,plot=F,out="geo.pdf") {
	suppressPackageStartupMessages(require(DESeq2))

	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}

	countData<- X$countData[rowSums(X$countData)>=cutoff,]
	dds <- 	DESeqDataSetFromMatrix(countData,X$colData,design)
	geoMeans = apply(counts(dds), 1, gm_mean)
	diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)
	if (plot) {
		desPCAplotter(dds,out)
	}
	
	return(DESeq(diagdds, fitType="local"))
}


desPCAplotter<- function(dds,out) {
	rld <- tryCatch( {
		print("Calculating VST")
		varianceStabilizingTransformation(dds,blind=F,fitType="local")
	}, error = function(e) {
		print("Unable to calculate VST")
		rld <- tryCatch( {
			DESeqTransform(SummarizedExperiment(log2(counts(estimateSizeFactors(dds), normalized=TRUE) + 1),colData=colData(dds)))
		}, error = function(e) {
			print("Unable to calculate DST")
			print("skipping rlog - PCA will use untransformed data") 
			return(dds)
		})
		return(rld)
	})
	print("plotting PCA")
	pdf(out,height=8,width=8)
	plotPCAWithLabels(rld)
	dev.off()
}

