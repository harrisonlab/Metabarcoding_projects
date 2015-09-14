library(biom)
library(sqldf)
library(DESeq2)

args <- commandArgs(TRUE) 
print("loading biom data\n")
bm.biom <- read_biom(args[1])

print("extracting useful data\n")
df.biom.data <- data.frame(as.matrix(biom_data(bm.biom)))
df.biom.taxon <- as.data.frame(do.call(rbind,observation_metadata(bm.biom, 1:length(bm.biom$rows))))
colnames(df.biom.taxon) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

df.biom <- merge(df.biom.taxon,df.biom.data,by="row.names",all.x=TRUE)
rownames(df.biom) <- df.biom$Row.names
df.biom<- df.biom[-1]

no_samples <- 7 + as.integer(args[2])
countData<- df.biom[(rowSums(df.biom[,8:no_samples])>=3),8:no_samples]
colData <- read.table("colData")

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
sig.res <- sqldf("select * from 'res.merge' where padj<=0.05",row.names=TRUE)
write.csv(sig.res,args[4])
