import_ubiom <- function (
	locX,
	locY,
	locZ
){
	options(stringsAsFactors = FALSE)
	countData <- read.table(locX,header=T,sep="\t", comment.char="")
	rownames(countData ) <- countData [,1]
	countData <- countData [,-1]
	taxonomy <- read.csv(locY,header=F)
	taxonomy <- taxonomy [,c(1,2,4,6,8,10,12,14)]
	rownames(taxonomy) <- taxonomy[,1]
	taxonomy <- taxonomy[,-1]
	colnames(taxonomy) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
	colData <- read.table(locZ,sep="\t",header=T)
	rownames(colData) <- colData [,1]
	colData <- colData[,-1,drop=FALSE]
	countData <- countData[,rownames(colData)]
	ls.biom <- list(countData,colData, taxonomy)
	names(ls.biom) <- c("countData","colData","taxonomy")
	return(ls.biom)
}

ubiom_to_des <- function(
	obj, 
	design=~1,
	fit=F,
	calcFactors=function(d)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(DESeq2))

	dds <- 	suppressWarnings(DESeqDataSetFromMatrix(obj[[1]],obj[[3]],design))

	sizeFactors(dds) <- calcFactors(dds)

    	if (fit) {
    	 	return(DESeq(dds,...))
    	} else {
    		return(dds)
    	}
} 	

ubiom_to_phylo <- function(obj){
	 phyloseq(
	 	otu_table(obj[[1]],taxa_are_rows=T),
	 	tax_table(as.matrix(obj[[2]])),
	 	sample_data(obj[[3]])
	 )
}

phylo_to_des <- function(
	obj,
	design=~1,
	fit=F,	
	calcFactors=function(d)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(DESeq2))
	dds <-  phyloseq_to_deseq2(obj,design)
	sizeFactors(dds) <- calcFactors(dds)
    	if (fit) {
    	 	return(DESeq(dds,...))
    	} else {
    		return(dds)
    	}
} 

phylo_to_ubiom <- function(obj) {
	list(
		countData=as.data.frame(obj@otu_table@.Data),
		taxonomy=as.data.frame(obj@tax_table@.Data),
		colData=as.data.frame(suppressWarnings(as.matrix(sample_data(obj)))) # suppresses a warning from the matrix call about loss of S4 state
	)
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


plotPCA <- function (	
	obj, 
	design = "condition",
	labelby,
	ntop = 500,
	pcx = 1,
	pcy = 2, 
	returnData = FALSE,
	cofix=F,
	trans=T,
	transform=function(o,design,...) 
	{	
		suppressPackageStartupMessages(require(DESeq2))
		dots <- list(...)
		if(!is.null(dots$calcFactors)) {
			calcFactors <- dots$calcFactors
			dots$calcFactors<-NULL
			if(length(dots)>1) {
				assay(varianceStabilizingTransformation(phylo_to_des(o,design,calcFactors=calcFactors),unlist(dots)))
			} else {
				assay(varianceStabilizingTransformation(phylo_to_des(o,design,calcFactors=calcFactors)))
			}
		} else {
			assay(varianceStabilizingTransformation(phylo_to_des(o,as.formula(design)),...))
		}
	},...
) {
	suppressPackageStartupMessages(require(genefilter))
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(DESeq2))

	if(trans) {
		obj@otu_table@.Data <- transform(obj,as.formula(paste0("~",design)),...)
	}
    
 	if (returnData) {
 		#d <- pca$x
 		#attr(d, "percentVar") <- percentVar
 		pca <- prcomp(t(otu_table(obj)))
 		pca$percentVar <- pca$sdev^2/sum(pca$sdev^2)
 		return(pca)
	}
 
 	rv <- rowVars(otu_table(obj))
 
	colData <- sample_data(obj)
	#suppressWarnings(as.data.frame(as.matrix(obj@sam_data)))
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	pca <- prcomp(t(otu_table(obj)[select, ]))
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	

	
	if (!all(design %in% names(colData))) {
		stop("the argument 'design' should specify columns of colData")
	}
	design.df <- as.data.frame(colData[, design,drop = FALSE])
	group <- if (length(design) > 1) {
		factor(apply(design.df, 1, paste, collapse = " : "))
	}
	else {
		as.factor(sample_data(obj)[[design]])
	}
	
	if (!missing(labelby)) {
		shape <- as.factor(sample_data(obj)[[labelby]])
		d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,design.df,shape = shape)
		colnames(d)[grep("shape", colnames(d))] <- labelby
	} else {
		d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,design.df)
	}

	colnames(d)[grep("group", colnames(d))] <- design

	if(cofix) {
		d[,1] <- d[,1] * percentVar[pcx]
		d[,2] <- d[,2] * percentVar[pcy]
	}

	g <- ggplot()
	g <- g + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
	if (!missing(labelby)) {
	g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group, shape=shape),size=3)
	g <- g + scale_shape_discrete(name=labelby)
	} else {
	g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=3)
	}
	g <- g + scale_colour_discrete(name=design)
	g <- g + xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance"))
	g <- g + ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
	return(g)
}



plotOrd <- function (	
	obj, 
	colData,
	design = "condition",
	shapes,
	labels=F,
	cluster=NULL,
	continuous=F, 
	xlims=NULL,
	ylims=NULL,
	legend=T,
	xlabel="Dimension 1",
	ylabel="Dimension 2",
	dimx=1,
	dimy=2,
	pointSize=2,
	...
) {

	suppressPackageStartupMessages(require(ggplot2))
   
	if (!all(design %in% names(colData))) {
		stop("the argument 'design' should specify columns of colData")
	}
	design.df <- as.data.frame(colData[, design,drop = FALSE])
	group <- if (length(design) > 1) {
		factor(apply(design.df, 1, paste, collapse = " : "))
	}
	else {
		if(continuous) {
			colData[[design]]
		} else {	
			as.factor(colData[[design]])
		}
	}
	
	if (!missing(shapes)) {
		shape <- as.factor(colData[[shapes]])
		d <- data.frame(PC1 = obj[, dimx], PC2 = obj[, dimy], group = group,design.df,shape = shape)
		colnames(d)[grep("shape", colnames(d))] <- shapes
	} else {
		d <- data.frame(PC1 = obj[, dimx], PC2 = obj[, dimy], group = group,design.df)
	}

	colnames(d)[grep("group", colnames(d))] <- design

	g <- ggplot()
	g <- g + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)
	g <- g + theme_bw()
	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
	if(!legend) {
		g <- g+ theme(legend.position="none")
	}	
	if (!missing(shapes)) {
		g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group, shape=shape),size=pointSize)
		g <- g + scale_shape_discrete(name=shapes)
	} else {
		g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=pointSize)
	}
	if(labels) {
		g <- g + geom_text(data=d, mapping=aes(x=PC1, y=PC2, label=row.names(obj),colour=group), size=(pointSize+1), vjust=2, hjust=0.5)
	}

	if(continuous) {
		#g <- g + scale_color_gradientn(colours = rainbow(5))
		g <- g + scale_color_gradient(low="red", high="green",name=design)
	} else {
		g <- g + scale_colour_discrete(name=design)

		if(!is.null(cluster)) {
			km <- kmeans(obj,...)
			d$Cluster<-km$cluster
			g<-g+stat_ellipse(data=d,mapping=aes(x=PC1,y=PC2,fill=factor(Cluster)), geom="polygon", level=cluster, alpha=0.2)
		}
	
	}
	g <- g + xlab(xlabel)
	g <- g + ylab(ylabel)
	return(g)
}



plotTaxa <- function(
	obj=mybiom, 	# obj (phloseq) a phyloseq object which must include taxonomy and sample data (or alternatively an S3 list)
	taxon="phylum", # taxon (str) is the taxonomic level of interest
	design, 	# condition (str) describes how the samples should be grouped (must be column of sample data)
	proportional=T,	# proportional (bool) whether the graph should use proportional or absolute values
	cutoff=1, 	# cutoff (double) for proportional graphs. 
	topn=0, 	# topn (int)taxons to display (by total reads) for non-prortional graphs. T
	others=T, 	# combine values less than cutoff/topn into group "other"
	ordered=F, 	# order by value (max to min)
	type=1, 	# type: (1) by sample stacked (2) by taxonomy stacked (3) by taxonomy compared
	fixed=F, 	# fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
	ncol=1, 	# ncol is a ggplot paramter to use n columns for the legend
	ret_data=F,
	coloured=T,
	trans=T,	# set to False if obj already contains pre-transformed counts (useful for larger OTU tables) 
	transform=function(o,design,...) 
	{	
		suppressPackageStartupMessages(require(DESeq2))
		dots <- list(...)
		
		if(!is.null(dots$calcFactors)) {
			calcFactors <- dots$calcFactors
			dots$calcFactors<-NULL
			if(length(dots)>=1) {
				assay(varianceStabilizingTransformation(ubiom_to_des(o,design,calcFactors=calcFactors),unlist(dots)))
			} else {
				assay(varianceStabilizingTransformation(ubiom_to_des(o,design,calcFactors=calcFactors)))
			}
		} else {
			assay(varianceStabilizingTransformation(ubiom_to_des(o,as.formula(design)),...))
		}
	}, # data transformation function 
	... # arguments to pass to transform function (obviously they could just be set in the function, but this looks neater)
) {
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(scales))
	
	if(isS4(obj)) {
		obj <- phylo_to_ubiom(obj)
	} 
	
	if(trans) {
		temp <- design
		idx <- grep(design,colnames(obj[[3]]))
		if(length(unique(obj[[3]][idx]))==1) {
			design<-1
		}
		obj[[1]] <- as.data.frame(transform(obj,as.formula(paste0("~",design)),...))
		design<-temp
	}
	#obj[[1]][obj[[1]]] <- obj[[1]][obj[[1]]]+abs(min(obj[[1]][obj[[1]]]))

	obj[[1]][obj[[1]]<0] <- 0
	
	taxa_sum <- sumTaxa(obj,taxon=taxon,design=design)
	taxa_sum$taxon[grep("\\(",taxa_sum$taxon)] <- taxa_sum$taxon[sub("\\(.*"," incertae sedis",taxa_sum$taxon)]

	if(!topn) {
		obj[[3]]$MLUflop <- 1 #assigns the MLU flop digit
		tx <- sumTaxa(obj,taxon=taxon,"MLUflop")
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

	taxa_cut <- taxa_sum[taxa_sum[,1]%in%txk,]
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=T),]
	if(others) {
		taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="others" ,t(colSums(taxa_sum[!taxa_sum[,1]%in%txk,-1]))),names(taxa_cut)))
	}
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	if(ordered) {
		taxa_cut[,ncol(taxa_cut)+1] <- 0
		taxa_cut[,1] <- reorder(taxa_cut[,1],-rowSums(taxa_cut[,-1]))
		taxa_cut <- taxa_cut[,-ncol(taxa_cut)]
	}
	if(ret_data) {
		return(taxa_cut)
	}
	md2 <- melt(taxa_cut,id=colnames(taxa_cut)[1])
	md2$variable <- factor(md2$variable, levels=levels(md2$variable)[order(levels(md2$variable))]  )
	if (type==1) {
		g <- ggplot(md2,aes_string(x=md2[,2],y=md2[,3],fill=taxon))
	} else if (type==2) {
		colnames(md2) <- c("taxa",design,"value")
		g <- ggplot(md2,aes_string(x=as.factor(md2[,1]),y=md2[,3],fill=design))
	}
	if(coloured) {
		g <- g + geom_bar(stat="identity",colour="white")

	 } else {
		g<-g+geom_bar(stat="identity",colour="white",fill="black")
	}
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


calcCorrelog<- function(pca,obj,pc,na.add,returnInp,returnCD) {
	cond<-"Y"

	pc.x <- scores(pca)[rownames(scores(pca))%in%rownames(sample_data(obj)[sample_data(obj)$condition==cond]),]
	col.x <- col.x <- sample_data(obj)[sample_data(obj)$condition==cond,]

	pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
	pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
	names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
	dvec <- pc.reshape$distance

	if (!missing(na.add)) {
		inp <- pc.reshape[[pc]]
		dvec<- unlist(sapply(1:length(inp),function(i) if(i%in%na.add){return(c(mean(c(dvec[i],dvec[(i-1)])),dvec[i]))}else{return(dvec[i])}))
		inp1 <- sapply(1:length(inp),function(i) if(i%in%na.add){return(c(NA,inp[i]))}else{return(inp[i])})
	}else {
		inp1<-pc.reshape[[pc]]
	}
	
	cond<-"N"

	pc.x <- scores(pca)[rownames(scores(pca))%in%rownames(sample_data(obj)[sample_data(obj)$condition==cond]),]
	col.x <- col.x <- sample_data(obj)[sample_data(obj)$condition==cond,]
	pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
	pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
	names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])	

	if (!missing(na.add)) {
		inp <- pc.reshape[[pc]]
		inp2 <- sapply(1:length(inp),function(i) if(i%in%na.add){return(c(NA,inp[i]))}else{return(inp[i])})
	}else {
		inp2<-pc.reshape[[pc]]
	}
	
	if(returnInp) {
		return(cbind(unlist(inp1),unlist(inp2),dvec))
	}

	ct1 <- correr1(unlist(inp1),returnCD)
	ca1 <- correr1(unlist(inp2),returnCD)

	d<-as.data.frame(cbind(ct1,ca1,dvec[1:(length(dvec)-2)]))
	d
}


plotCorrelog <- function(
	pca,
	obj,
	pc="PC1",
	cutoff=15,
	xlim=NULL,
	ylim=NULL,
	na.add,
	returnData=F,
	returnInp=F,
	returnCD=F,
	data,
	legend=T,
	lineWidth=1.5 
) {

	if(returnInp|returnCD){
		returnData=T
	}
	
	if(!missing(data)){
		d<-data
	} else {
		d <- calcCorrelog(pca,obj,pc,na.add,returnInp,returnCD)
	}
	
	if(returnData) {
		return(d)
	}
	
	d<- d[1:(length(d$V3[d$V3<=cutoff])),]
	names(d) <- c(paste0("Tree_",pc),paste0("Aisle_",pc),"Distance")
	d2 <- melt(d,id="Distance")
	colnames(d2)[3] <- c("Correlation")

	g <- ggplot(d2)
	g <- g + coord_fixed(ratio = 10, xlim = xlim, ylim = ylim, expand = TRUE)
	g <- g + theme_bw()
	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
	if(!legend) {
		g <- g+ theme(legend.position="hidden")
	}
	g <- g + geom_line(na.rm=T,aes(x=Distance, y=Correlation, colour=variable),size=lineWidth)
	g <- g + scale_colour_manual(values=c("red","green"))
	g

}



phyloTaxaTidy <- function(obj) {
	colnames(obj) <- c("kingdom","phylum","class","order","family","genus","species")
	obj <- sub("*._+","",obj)
	obj <- t(apply(obj,1,taxonomyTidy))
	return(obj)
}

taxonomyTidy <- function(x) {
	if (x[2]=="unknown") {x[2] <- paste(x[1],"(k)",sep="")}
	if (x[3]=="unknown") {if(any(grep('\\(',x[2]))) {x[3]<-x[2]}else{x[3]<-paste(x[2],"(p)",sep="")}}
	if (x[4]=="unknown") {if(any(grep('\\(',x[3]))) {x[4]<-x[3]}else{x[4]<-paste(x[3],"(c)",sep="")}}
	if (x[5]=="unknown") {if(any(grep('\\(',x[4]))) {x[5]<-x[4]}else{x[5]<-paste(x[4],"(o)",sep="")}}
	if (x[6]=="unknown") {if(any(grep('\\(',x[5]))) {x[6]<-x[5]}else{x[6]<-paste(x[5],"(f)",sep="")}}
	if (x[7]=="unknown") {if(any(grep('\\(',x[6]))) {x[7]<-x[6]}else{x[7]<-paste(x[6],"(g)",sep="")}}
	return(x)
}

sumTaxa <- function(
	obj,
	taxon="phylum",
	design="condition",
	proportional=F
){
# sums by sample data 
	suppressPackageStartupMessages(require(plyr))
	suppressPackageStartupMessages(require(reshape2))
	tx <- obj[[2]][,taxon]
	dtx <- cbind(obj[[1]],tx)
	md <- melt(dtx,id="tx")
	obj[[3]]$all <- "all"
	md$variable <- mapvalues(md$variable,from=rownames(obj[[3]]), to=as.character(obj[[3]][,design]))
	nd <- dcast(md,...~variable,sum)
	colnames(nd)[1] <- taxon
	if(proportional) {
		nd[-1] <-  prop.table(as.matrix(nd[,-1]),2)*100
	}
	return(nd)
}

topTaxa <- function(
	obj
)


combine_biom <- function(locX,locY) {
	biom1 <- read.table(locX,header=T,sep="\t", comment.char="")	
	biom2 <- read.table(locY,header=T,sep="\t", comment.char="")
	biom <- merge(biom1,biom2,by.x="X.OTU.ID",by.y="X.OTU.ID",all=T)
	biom[is.na(biom)] <- 0
	return(biom)	
}

fltTaxon <- function(
	obj,
	taxon="phylum",
	out="phylo"
){
# same as phyloseq tax_glom (drops NA columns returned by tax_glom), but works on S3 biom data (i.e. ubiom)
# perhaps tax_glom is good with big datasets as fltTaxon is miles faster - aggregate will get pretty slow for large datasets
	if(class(obj)[[1]]=="phyloseq") {
		obj <- phylo_to_ubiom(obj)
	}
	n <- which(colnames(obj[[2]])==taxon)
	x <- aggregate(obj[[1]],by=obj[[2]][,1:n],sum)
	ls.biom <- list(x[,(n+1):ncol(x)],x[,1:n],obj[[3]])
	names(ls.biom) <- c("countData","taxonomy","colData")
	if(out=="phylo") {
		return(ubiom_to_phylo(ls.biom))
	}
	return(ls.biom)	
}

# function sometimes useful for replacing calcfactors
geoMeans <- function(d) {
	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}	
	gm = apply(counts(d), 1, gm_mean)
	sizeFactors(estimateSizeFactors(d, geoMeans =gm))
}

correr1 <- function(x,returnData=F) {
	y <- x
	count<-1
	mycorr <- NULL
	while (length(x) >2) {
		if(returnData) {
			mycorr[[count]] <- cbind(y[1:(length(x)-1)],x[2:length(x)])
		} else {
			mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)],use="pairwise.complete.obs")
		}
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}



correr2 <- function(X,breaks) {
	y <- seq(0,max(X[,2])+1,breaks)
	vec <- X[,2] 
	test <- nearest.vec(y,vec)	
	test[9] <- NA
	x <- X[X[,2]==test,1]
	mycorr <- numeric(0)
	count<-1
	while (length(x) >2) {
		mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)],use="pairwise.complete.obs")
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}

nearest.vec <- function(x, vec,breaks=1) {
    smallCandidate <- findInterval(x, vec, all.inside=TRUE)
    largeCandidate <- smallCandidate + 1
    nudge <- 2 * x > vec[smallCandidate] + vec[largeCandidate]
    empty <- abs(vec[smallCanditate]-break)
    return(smallCandidate + nudge)
}



correr1_dat <- function(x) {
	y <- x
	mycorr <- NULL
	count<-1
	while (length(x) >2) {
		mycorr[[count]] <- cbind(y[1:(length(x)-1)],x[2:length(x)])
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}