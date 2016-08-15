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

Set strinAsFactor to false - apparently this can be set in .Rprofile, but this doesn't work on my production cluster environment
```{r}
options(stringsAsFactors = FALSE)
```
### phyloseq
Import biom table, taxonomy and sample data into a phyloseq object
Create phyloseq object in R
```{r}
library("phyloseq")
biom_file


mybiom <- import_biom("16S.taxa.biom",refseqfilename="16S.otus.fa")
```


### DESeq2


Requires analysis2.R and deseq.r

ubiom makes a S3 biom object from the OTU table (16S.otu_table.txt), OTU taxonomy (16S.taxa) and sample description file (colData) analysis2.R/deseq.r contain scripts to produce deseq objects and run differential analysis + a few graphing options.

The OTU table header starts with a hash. To import into R set comment.char="" in the read.table parameter
