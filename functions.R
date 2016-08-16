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

