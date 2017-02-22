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
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c NEMpre \
 "$ARDERI/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $ARDERI/data/$RUN/$SSU/fasta \
 $ARDERI/metabarcoding_pipeline/primers/nematode.db \
 $MINL $MAXL $QUAL
```

#### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2

```shell
mkdir -p $ARDERI/data/$RUN/$SSU/filtered
find $ARDERI/data/$RUN/$SSU/fasta -type f -name *.r*|xargs -I myfile mv myfile $ARDERI/data/$RUN/$SSU/filtered/.

cd $ARDERI/data/$RUN/$SSU/filtered
for f in $ARDERI/data/$RUN/$SSU/filtered/*r1.fa
do
    R1=$f
    R2=$(echo $R1|sed 's/\.r1\.fa/\.r2\.fa/')
    S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
    $ARDERI/metabarcoding_pipeline/scripts/catfiles_v2.pl $R1 $R2 $S;
done

mkdir R1
mkdir R2
mv *r1* R1/.
mv *r2* R2/.
```

## UPARSE

### Cluster 

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE \ $ARDERI $RUN $SSU 0 0
```
### Assign taxonomy

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign \ $ARDERI $RUN $SSU 
```

### Create OTU tables

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU \ $ARDERI $RUN $SSU $FPL $RPL true
```


###[16S workflow](../master/16S%20%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)
