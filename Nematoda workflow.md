# Nematode workflow


## Pre-processing
Script will:<br>
1. Remove reads with both forward and reverse primers<br>
2. Remove reads with adapter contamination<br>
3. Filter for quality and minimum length (with UTRIM)<br>
4. Convert FASTQ to single line FASTA

```shell

$METAGENOMICS/scripts/ARDERI.sh -c ITSpre \ 
 $METAGENOMICS/data/$RUN/ITS/fastq/*R1*.fastq \ 
 $METAGENOMICS/data/$RUN/ITS/fasta \
 $METAGENOMICS/primers/primers.db \
 200 200 1; 
```

### SSU/58S/LSU removal 

#### Identify SSU, 5.8S  and LSU regions

This will create a large number of array jobs on the cluster

Fungi
```shell
$METAGENOMICS/scripts/ARDERI.sh -c procends \
 $METAGENOMICS/data/$RUN/$SSU/fasta \
 R1 \
 $METAGENOMICS/hmm/ssu_end.hmm 	\
 $METAGENOMICS/hmm/58s_start.hmm \
 ssu 58ss 20

$METAGENOMICS/scripts/ARDERI.sh -c procends \
 $METAGENOMICS/data/$RUN/$SSU/fasta \
 R2 \
 $METAGENOMICS/hmm/lsu_start.hmm \
 $METAGENOMICS/hmm/58s_end.hmm \
 lsu 58se 20


```

Oomycetes
```shell
$METAGENOMICS/scripts/ARDERI.sh -c procends \
 $METAGENOMICS/data/$RUN/$SSU/filtered \
 "" \
 $METAGENOMICS/hmm/others/Oomycota/ssu_end.hmm \
 $METAGENOMICS/hmm/others/Oomycota/58s_start.hmm \
 ssu 58ss 20
```

#### Remove SSU, 5.8S  and LSU regions and merge output

If reverse read quality was poor and it was necessary to truncate reads to get more than a couple of reads past set LOWQUAL to TRUE

LOWQUAL keeps reads which lack 5.8S homology - this is necessary as trimming will in most instances have removed the homologous region

Fungi
```shell
$METAGENOMICS/scripts/ARDERI.sh -c ITS \
 "$METAGENOMICS/data/$RUN/$SSU/fasta/*R1" \
 $METAGENOMICS/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58"

LOWQUAL=FALSE   
$METAGENOMICS/scripts/ARDERI.sh -c ITS \
 "$METAGENOMICS/data/$RUN/$SSU/fasta/*R2" \
 $METAGENOMICS/scripts/rm_58Se_LSU_v2.R \
 "*.\\.58" \
 "*.\\.lsu" \
 $LOWQUAL
```

Oomycetes
```shell
$METAGENOMICS/scripts/ARDERI.sh -c ITS "$METAGENOMICS/data/$RUN/$SSU/filtered/*D" $METAGENOMICS/scripts/rm_SSU_58Ss.R "*.\\.ssu" "*.\\.58"
```

There's a slight problem with one of the scripts and the fasta names...
```shell
for f in *.fa; do
	sed -i -e 's/ .*//' $f
done
```



#### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2

```shell
mkdir -p $METAGENOMICS/data/$RUN/$SSU/filtered
find $METAGENOMICS/data/$RUN/$SSU/fasta -type f -name *.r*|xargs -I myfile mv myfile $METAGENOMICS/data/$RUN/$SSU/filtered/.

cd $METAGENOMICS/data/$RUN/$SSU/filtered
for f in $METAGENOMICS/data/$RUN/$SSU/filtered/*r1.fa
do
    R1=$f
    R2=$(echo $R1|sed 's/\.r1\.fa/\.r2\.fa/')
    S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
    $METAGENOMICS/scripts/catfiles_v2.pl $R1 $R2 $S;
done

mkdir R1
mkdir R2
mv *r1* R1/.
mv *r2* R2/.
```

## UPARSE
FPL=23 
RPL=21

### Cluster 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

```shell
$METAGENOMICS/scripts/ARDERI.sh -c UPARSE \ $METAGENOMICS $RUN $SSU 0 0
```
### Assign taxonomy
NOTE:- I still need to build nematode utax taxonomy database from Silva_SSU.

```shell
$METAGENOMICS/scripts/ARDERI.sh -c tax_assign \ $METAGENOMICS $RUN $SSU 
```

### Create OTU tables

Concatenates unfiltered reads, then assigns forward reads to OTUs. For any non-hits, attemps to assign reverse read (ITS2) to an OTU. 

```shell
$METAGENOMICS/scripts/ARDERI.sh -c OTU \ $METAGENOMICS $RUN $SSU $FPL $RPL true
```


###[16S workflow](../master/16S%20%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)
