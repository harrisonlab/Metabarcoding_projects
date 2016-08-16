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
```{r}

# res <- estimate_richness(mybiom) ## data used for plot_richness graphs
pdf("16S.alpha_bysex.pdf", height=8,width=8)
plot_richness(mybiom,x="condition",color="Sex",measures=c("Chao1", "ACE", "Shannon", "Simpson"))
dev.off()
```
##### beta diversity
```{r}
pdf("16S.beta-diversity.pdf",height=8,width=8)
plotPCA(dds)
dev.off()
```
#### taxa graphs
```{r}
pdf("16S.phylum.pdf",height=8,width=8)
plotTaxa(mybiom,"phylum","condition")
dev.off()
```

