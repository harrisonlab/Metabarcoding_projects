#Statistical analysis

## Combine biom and taxa
biom_make.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biom and output a combined taxa and biom file to standard out

e.g. for 16S
```shell
cat 16S.rdp|$METAGENOMICS/scripts/mod_taxa.pl >16S.taxa
$METAGENOMICS/scripts/biom_maker.pl 16S.taxa 16S.otu_table.biom >16S.taxa.biom
```

## R analysis
R will import table with character columns as factors. This might be prefereable if your doing lots of linear modeling, but for any data manipulation it is a disaster and will lead to unexpected (or worse unnoticed) errors.

Set stringAsFactor to false - apparently this can be set in .Rprofile, but this doesn't work on my production cluster environment
```{r}
options(stringsAsFactors = FALSE)
```
There are several R functions throughout this module which are found in functions.R
```{r}
source("../../scripts/functions.R")
```

### phyloseq
Import biom table, taxonomy and sample data into a phyloseq object
Create phyloseq object in R

The taxonomy imported from UPARSE doesn't include rank names and adds extra stuff to the taxa names (k__,p__ and etc.)
phyloTaxaTidy will fix this and change "unknown" to the lowest known rank (and append a character indicating rank)
```{r}
library("phyloseq")
biom_file = "16S.taxa.biom"
otu_file = "16S.otus.fa" # might be useful at some stage
colData = "colData"
mybiom <- import_biom(biom_file) # ,refseqfilename=out_file
sample_data(mybiom) <- read.table(colData,header=T,sep="\t",row.names=1)
tax_table(mybiom) <- phyloTaxaTidy(tax_table(mybiom))

# An example of removing certain OTUs from a phyloseq object - and using S4 methods to access the data.
# This will filter based on OTU present in the first column "condition" of colData 
t1 <- aggregate(t(otu_table(mybiom)),by=list(sample_data(mybiom)[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>0,mybiom)
```

##### Merge ITS1 and ITS2
```{r}
ITS1 = "ITS1.taxa.biom"
ITS2 = "ITS2.taxa.biom"
mybiom <- merge_phyloseq(import_biom(ITS1),import_biom(ITS2))
```

##### Create and add phylogentic tree to mybiom

```{r}
library(ape)
temp_connection = file("ITS.phy", 'r')
len = readLines(temp_connection, n=1)
len = as.numeric(len)
close(temp_connection)
phylip_data = read.table("ITS.phy", fill=T, row.names=1, skip=1, col.names=1:len)
ITS.nj <- nj(as.dist(phylip_data))
write.tree(ITS.nj,"ITS.tree")
phy_tree(mybiom) <- ITS.nj
```

#### Beta-diversity statistical analysis
Using PERMANOVA (This needs additional work, can't find any info on whether to use normised or raw reads - if bray-curtis is non-parametric then no need for library size normalisation, but I'm using scale to transform by translation and expansion what's the point of doing this??)
```{r}
library(vegan)
obj <- mybiom
obj@otu_table@.Data <- counts(phylo_to_des(mybiom,design=~1),normalize=T)
d <- as.data.frame(cbind(sample_data(obj),t(otu_table(obj))))
euclid <- scale(d[,c(-1,-2)]) # minus however many columns colData has defined
adonis(euclid~condition,d,method='bray')
```

### Spatial analysis
Simple first step - correspondence analysis

1. All data
2. Tree station
3. Aisle

```{r}
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]=="experiment",mybiom)

f1 <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
f2 <- prune_samples(sample_data(myfiltbiom)[[1]]=="Y",myfiltbiom)
f3 <- prune_samples(sample_data(myfiltbiom)[[1]]=="N",myfiltbiom)

f1 <- prune_taxa(rowSums(otu_table(f1))>2,f1)
f2 <- prune_taxa(rowSums(otu_table(f2))>2,f2)
f3 <- prune_taxa(rowSums(otu_table(f3))>2,f3)

CCA1 <- ordinate(f1,method="CCA",formula=~distance)
CCA2 <- ordinate(f2,method="CCA",formula=~distance)
CCA3 <- ordinate(f3,method="CCA",formula=~distance)

plot(CCA1)
plot(CCA2)
plot(CCA3)
dev.off()

anova(CCA1) # permutation analysis
anova(CCA2)
anova(CCA3)

```
More advanced technique; Principal Coordinates of Neighbour Matrices
```{R}
library(vegan)

# apply various filters to data
myfiltbiom <- prune_samples((sample_data(mybiom)[[10]]=="experiment")&(sample_data(mybiom)[[1]]!="C"),mybiom)
myfiltbiom@sam_data$gap[myfiltbiom@sam_data$condition=="N"] <- 0
myfiltbiom <- prune_taxa(rowSums(otu_table(myfiltbiom))>2,myfiltbiom)

#calculate euclidean distance between sample points
euclid <- dist(sample_data(myfiltbiom)[,6:7],method="euclidean")
#unweighted PCNM
#pcnm1 <- pcnm(euclid)
#calculate PCNM weights
cs <- colSums(otu_table(myfiltbiom))/sum(otu_table(myfiltbiom))
#weighted PCNM (this is more useful for cca)
pcnm2 <- pcnm(euclid,w=cs)
#cca on OTU data with positive eigen vectors from PCNM as the independent variables
#the residual should have no distance trend (i.e. the CA component of the ord object)
ord <- cca(t(otu_table(myfiltbiom))~scores(pcnm2))
anova(ord)
plot(ord)
msoplot(mso(ord, sample_data(myfiltbiom)[,6:7]))
dev.off()
#spatial free data
newCounts <- t(ord$CA$Xbar)
## 
ord2 <- cca(t(otu_table(myfiltbiom))~scores(pcnm2)[,1]+scores(pcnm2)[,2]+scores(pcnm2)[,3..n])
anova(ord2)
anova(ord2,by="terms",parellel=12,model="direct",permu=2000)
```


```{R}
obj <- mybiom
obj@otu_table@.Data <- assay(varianceStabilizingTransformation(phylo_to_des(obj)))
mynmds <- ordinate(obj,method = "NMDS",distance="bray",autotransform=F,try=100)
#plotOrd(mynmds$points,sample_data(obj),design="condition",shapes="location")
plot_ordination(obj,mynmds,color="condition",shape="location")
```



### DESeq2
It's possible to convert a phyloseq object to a DESeq datamatrix with the wrapper function phylo_to_des.phylo_to_des has the option to specify the size factor calculation using the option calcFactors (see plotTaxa for examples of how to use this option). Set fit=T to fit a GLM model to the data. Further arguments will be passed to DESeq.
```{r}
dds <- phylo_to_des(mybiom,fit=T, fitType="local",...)
```
#### Differential OTU abundance
Using DESeq2 it's possible to calculate the probability of OTUs having different abundances between condtions. The default will use the the condition column of the dds object's colData table, and take the first two conditions. To specify a different column or use different "condtions use the contrast=c("column_name","condition_1","condition_2") construct when calling the results method.

This code will merge the phyloseq taxonomy object to the results.

```{r}
alpha <- 0.05 # significance level
res = results(dds, alpha=alpha)	
#res = results(dds, alpha=alpha,contrast=c("condition","N","K"))	## specify different contrasts to make
res.merge <- merge(as.data.frame(res),tax_table(mybiom),by="row.names",all.x=TRUE)
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
NMDS etc. plots 

```{R}
obj <- mybiom
obj@otu_table@.Data <- assay(varianceStabilizingTransformation(phylo_to_des(obj)))
mynmds <- ordinate(obj,method = "NMDS",distance="bray",autotransform=F,try=100)
#plotOrd(mynmds$points,sample_data(obj),design="condition",shapes="location")
plot_ordination(obj,mynmds,color="condition",shape="location")
```

plotOrd is a ggplot wrapper that does something similar to the phyloseq plotting method, but without the background and with axes on the same scale. Needs modifying to give a bit more control.



plotPCA is modified version of the DESeq2 version. 
It take the following options:

1. object (DESeq2 - required) a DESeq object 
2. intgroup (string - optional, default="condition") a column of colData used to describe (colour) the samples (e.g. infected/control)
3. labelby (string - optional) a 2nd column of colData used to descibe (shape) the samples (e.g. male/female)
4. ntop (int - optional, default=500) number of OTUs in descending count order to use in the PCA calculation
5. pcx (int - optional, default=1) X-axis pricipal component
6. pcy (int - optional, default=2) Y-axis pricipal component
7. returnData (bool - optional, default=F) returns the PCA results and exits
8. cofix (bool - optional, default=F) produces a graph with axes on the same scale
9. transform(fun - optional) a user supplied function to replace DESeq2 variance stabilising transform to transform the count matrix. A DESeq2 object will be passed to this function. 

plotPCAWithLabels will produce the same graph but with the addition of sample labels - useful for getting the name of outliers
```{r}
pdf("16S.beta-diversity.pdf",height=8,width=8)
plotPCA(dds)
dev.off()
```


#### taxa graphs
plotTaxa produces a ggplot2 bar chart of taxa counts
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
12. transform (fun - optional) a user supplied function to replace DESeq2 variance stabilising transform to transform the count matrix. An S3 list of countData, Taxonomy and colData is passed to the transform function . 
13. ... arguments to pass to transform function (o.k. these could just be set in the function, but this is a neater solution if the default function is used)

The default transform firts creates a DESeq2 object using ubiom_to_des. If all OTUs contain at least one zero value DDS can't calculate sizeFactors. The transform function could be modified to calculate VST on a data matrix, or a function "calcFactors" can be passed to the ubiom_to_des method which can calculate the sizeFatcors (e.g, edgeR calcNormFactors) 

```{r}
pdf("16S.phylum.pdf",height=8,width=8)
plotTaxa(mybiom,"phylum","condition")
plotTaxa(mybiom,"phylum","condition",blind=F) # passes  blind=F to transform function

# all OTUs contain at least one zer - pass calcFactors to ubiom_to_des
t1 <- aggregate(t(otu_table(mybiom)),by=list(sample_data(mybiom)[[1]]),FUN=sum)[-1]
myfiltbiom <- prune_taxa(apply(t1,2,prod)>==0,mybiom) # Keep OTUs with at least one zero value 
plotTaxa(myfiltbiom,design="condition",calcFactors=function(d){library(edgeR);sizeFactors(d) <- calcNormFactors(counts(d))}) 

# precalculate transform
dds <- phylo_to_des(mybiom,design=~1)
rld <- varianceStabilizingTransformation(dds)
mybiom <- mytransbiom
mytransbiom@otu_table@.Data <- assay(rld) # couldn't get method otu_table(mybiom) to accept rld data..
plotTaxa(mytransbiom,taxon="genus",design="condition",trans=F)
dev.off()
```

###[16S workflow](../master/16S%20%20workflow.md)
###[ITS workflow](../master//ITS%20workflow.md)
