# Nematode workflow

## Conditions
SSU determines the file location
FPL is forward primer length
RPL is reverse primer length

```shell
# nematodes
SSU=NEM 
FPL=23
RPL=18

# all
MINL=200
MAXL=300
QUAL=1

```

## Pre-processing
Script will:<br>
1. Remove reads with both forward and reverse primers<br>
2. Remove reads with adapter contamination<br>
3. Filter for quality and minimum length (with UTRIM)<br>
4. Convert FASTQ to single line FASTA

```shell
$METAGENOMICS/scripts/ARDERI.sh -c ITSpre \
 "$METAGENOMICS/data/$RUN/$SSU/fastq/*R1*.fastq" $RUN \
 $METAGENOMICS/data/$RUN/$SSU/fasta \
 $METAGENOMICS/primers/primers_nem.db \
 $MINL $MAXL $QUAL
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

### Cluster 

```shell
$METAGENOMICS/scripts/ARDERI.sh -c UPARSE \ $METAGENOMICS $RUN $SSU 0 0
```
### Assign taxonomy

```shell
$METAGENOMICS/scripts/ARDERI.sh -c tax_assign \ $METAGENOMICS $RUN $SSU 
```

### Create OTU tables

```shell
$METAGENOMICS/scripts/ARDERI.sh -c OTU \ $METAGENOMICS $RUN $SSU $FPL $RPL true
```


###[16S workflow](../master/16S%20%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)
