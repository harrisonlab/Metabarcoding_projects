# Oomycota workflow

## Conditions
SSU determines the file location
FPL is forward primer length
RPL is reverse primer length

```shell
# oomycota
SSU=OO 
FPL=21
RPL=20

# all
MINL=300
MINOVER=5
QUAL=0.5
```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size min_join_overlap max_errrors 

```shell
$METAGENOMICS/scripts/ARDERI.sh -c 16Spre \
	$METAGENOMICS/data/$RUN/$SSU/fastq/*R1*.fastq \
	$METAGENOMICS/data/$RUN/$SSU/filtered \
	$METAGENOMICS/primers/adapters.db \
	$MINL $MINOVER $QUAL
```
### SSU and 5.8S removal 

#### Identify SSU and 5.8S regions

This will create a large number of array jobs on the cluster

```shell
$METAGENOMICS/scripts/ARDERI.sh -c procends \
 $METAGENOMICS/data/$RUN/$SSU/filtered \
 "" \
 $METAGENOMICS/hmm/others/Oomycota/ssu_end.hmm \
 $METAGENOMICS/hmm/others/Oomycota/58s_start.hmm \
 ssu 58ss 20
```

#### Remove identified SSU and 5.8S regions

```shell
$METAGENOMICS/scripts/ARDERI.sh -c ITS \
	"$METAGENOMICS/data/$RUN/$SSU/filtered/*D" \
	$METAGENOMICS/scripts/rm_SSU_58Ss.R \
	"*.\\.ssu" "*.\\.58"
```

There's a slight problem with one of the scripts and the fasta names...
```shell
for f in *.fa; do
	sed -i -e 's/ .*//' $f
done
```

## OTU assignment 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

### Cluster 
```shell
$METAGENOMICS/scripts/ARDERI.sh -c UPARSE $METAGENOMICS $RUN $SSU 0 0
```
### Assign taxonomy
```shell
$METAGENOMICS/scripts/ARDERI.sh -c tax_assign $METAGENOMICS $RUN $SSU 
cp $SSU.otus.fa $SSU_v2.otus.fa
$METAGENOMICS/scripts/ARDERI.sh -c tax_assign $METAGENOMICS $RUN $SSU_v2
rm $SSU_v2.otus.fa
```

### Create OTU tables
```shell
$METAGENOMICS/scripts/ARDERI.sh -c OTU $METAGENOMICS $RUN $SSU $FPL $RPL
```


###[16S workflow](../master/16S%20%20workflow.md)
###[ITS workflow](../master//ITS%20workflow.md)
###[Nematoda workflow](../master/Nematoda%20%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)



