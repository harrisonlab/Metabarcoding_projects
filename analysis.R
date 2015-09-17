##############
#
# R script to analyses Quiim Biom data	
# Outputs: number of OTUS per sample, number of sequences per OTU 
# the number of unique taxa the OTUs are assigned to  per sample
# and number of statistically significant differences in OTU between sample groups
#
# Usage:
#
# Rscript analysis.R biom_file experimental_details median/geomeans output_file 	
#
############

library(biom)
library(DESeq2)

args <- commandArgs(TRUE) 

###biom
print("loading biom data\n")
bm.biom <- read_biom(args[1])
print("extracting useful data\n")
df.biom.data <- data.frame(as.matrix(biom_data(bm.biom)))
df.biom.taxon <- as.data.frame(do.call(rbind,observation_metadata(bm.biom, 1:length(bm.biom$rows))))
colnames(df.biom.taxon) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
df.biom <- merge(df.biom.taxon,df.biom.data,by="row.names",all.x=TRUE)
rownames(df.biom) <- df.biom$Row.names
df.biom<- df.biom[-1]

###DESSeq
countData<- df.biom[(rowSums(df.biom[,8:ncol(df.biom)])>=3),8:ncol(df.biom)]
colData <- read.table(args[2])
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)

if (args[3]=="median") {
	dds <- DESeq(dds, fitType="local")
	res = results(dds)
} else {

	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}

	geoMeans = apply(counts(dds), 1, gm_mean)
	diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)
	dds = DESeq(diagdds, fitType="local")
	res = results(dds)
}

res.merge <- merge(as.data.frame(res),df.biom,by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=0.05)
write.csv(sig.res,args[4])

print("no. of OTUs")
print(apply(df.biom.data,2,function(x) length(x[x>0])))
print("Sequences per OTU")
print(colSums(df.biom.data))
print("No. taxons")
X <- aggregate(df.biom[,8:ncol(df.biom)],by=list(df.biom$kingdom,df.biom$phylum,df.biom$class,df.biom$order,df.biom$family,df.biom$genus,df.biom$species),FUN=sum)
print(apply(X[,8:ncol(df.biom)],2,function(x) length(x[x>0])))
