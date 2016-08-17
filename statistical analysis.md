#Statistical analysis

## Combine biom and taxa
biom_make.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biom and output a combined taxa and biom file to standard out

e.g. for 16S
```shell
cat 16S.rdp|./mod_taxa.pl >16S.taxa
./biom_maker.pl 16S.taxa 16S.otu_table.biom >16S.taxa.biom
```

## R analysis
R will import table with character columns as factors. This might be prefereable if your doing lots of linear modeling, but for any data manipulation it is a disaster and will lead to unexpected (or worse unnoticed) errors.

Set stringAsFactor to false - apparently this can be set in .Rprofile, but this doesn't work on my production cluster environment
```{r}
options(stringsAsFactors = FALSE)
```
There are several R functions throughout this module which are found in functions.R
```{r}
source("functions.R")
```

### phyloseq
Import biom table, taxonomy and sample data into a phyloseq object
Create phyloseq object in R
```{r}
library("phyloseq")
biom_file = "16S.taxa.biom"
otu_file = "16S.otus.fa" # might be useful at some stage
colData = "colData"
mybiom <- import_biom(biom_file,refseqfilename=out_file)
sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)

# an example of removing certain OTUs from a phyloseq object
# this will filter based on OTU present in all conditions
t1 <- aggregate(t(mybiom@otu_table),by=list(mybiom@sam_data$condition),FUN=sum)
t1 <- t1[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>0,mybiom)

```
### DESeq2
It's possible to convert a phyloseq object to a DESeq datamatrix with phyloseq_to_deseq2
```{r}
dds <- phylo_to_des(mybiom)
```
#### Differential OTU abundance
Using DESeq2 it's possible to calculate the probability of OTUs having different abundances between condtions. The default will use the the condition column of the dds object's colData table, and take the first two conditions. To specify a different column or use different "condtions use the contrast=c("column_name","condition_1","condition_2") construct when calling the results method.

This code will merge the phyloseq taxonomy object to the results.

```{r}
alpha <- 0.05 # significance level
res = results(dds, alpha=alpha)	
#res = results(dds, alpha=alpha,contrast=c("condition","N","K"))	## specify different contrasts to make
taxa=mybiom@tax_table
colnames(taxa) <- c("kingdom","phylum","class","order","family","genus","species")
taxa <- sub("*._+","",taxa)
res.merge <- merge(as.data.frame(res),taxa,by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)
```
### plots

#### alpha diversity
Phyloseq has a builtin method for plotting alpha diversity graphs (plot_richness). Estimate_richness will return the data points used to plot the graphs (I may edit this function as the graphs are pretty basic - ugly anyway)

```{r}
# res <- estimate_richness(mybiom) ## data used for plot_richness graphs
pdf("16S.alpha_bysex.pdf", height=8,width=8)
plot_richness(mybiom,x="condition",color="Sex",measures=c("Chao1", "ACE", "Shannon", "Simpson"))
dev.off()
```
#### beta diversity
Beta diversity is plotted with a modified version of the DESeq2 plotPCA method. 
It take the following options:

1. object (DESeq2 - required) a DESeq object with size factors
2. intgroup (string - required, default="condition") a column of colData used to describe (colour) the samples (e.g. infected/control)
3. labelby (string - optional) a 2nd column of colData used to descibe (shape) the samples (e.g. male/female)
4. ntop (int - required, default=500) number of OTUs in descending count order to use in the PCA calculation
5. pcx (int - required, default=1) X-axis pricipal component
6. pcy (int - required, default=2) Y-axis pricipal component
7. returnData (bool - optional, default=F) returns the PCA results and exits
8. cofix (bool - optional, default=F) produces a graph with axes on the same scale
9. transform(fun - required, default VST) a function which describes how to transform the DDS size factors for plotting. Object will be passed to this function as its first option. 

plotPCAWithLabels will produce the same graph but with the addition of sample labels - useful for getting the name of outliers
```{r}
pdf("16S.beta-diversity.pdf",height=8,width=8)
plotPCA(dds)
dev.off()
```
I may update this to accept a phlyoseq object rather than DESeq and do the size factor calculations internally
#### taxa graphs
Produces a ggplot2 bar chart of taxa counts
It takes the following options:

1. obj (phyloseq - required) must is a phyloseq object which must include taxonomy and sample data
2. taxon (str - required) is the taxonomic level of interest
3. condition (str - required) describes how the samples should be grouped (must be column of sample data)
4. proportional (bool - optional, def=T) whether the graph should use proportional or absolute values
5. cutoff (double - optional, def =1) for proportional graphs. 
6. topn (int - optional)taxons to display (by total reads) for non-prortional graphs. 
7. others (bool - optional, def=T), combine values less than cutoff/topn into group "other"
8. reorder (bool - optional, =F) order by value (max to min)
9. type (int (1/2) - required, def=1) Type 1 produces a stacked barchart by sample, type 2 a barchart by taxa 
10. fixed (bool - optional, def=F) fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
11. ncol (int - optional, def=1) ncol is a ggplot paramter to use n columns for the legend
12. calcFactors (fun - optional) user supplied function to replace DESeq2 estimateSizeFactors, to calculate size factors on the OTU count matrix (DESeq2 object will be passed to function). estimateSizeFactors can't handle 0 count values for all samples. Something like the below will fix this error, or use a function (edgeR - calcNormFactors) which doesn't mind 0 counts.
function(d) {
	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}
	geoMeans = apply(counts(d), 1, gm_mean)
	sizeFactors(estimateSizeFactors(d, geoMeans = geoMeans))
}
13. transform (fun - optional) a function which describes how to transform the DESeq2 size factors for plotting (default is VST). DESeq2 object will be passed to this function.

```{r}
pdf("16S.phylum.pdf",height=8,width=8)
plotTaxa(mybiom,"phylum","condition")
dev.off()
```
