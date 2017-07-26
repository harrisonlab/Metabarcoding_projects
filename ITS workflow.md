# ITS workflow
```shell
SSU=ITS
FPL=23 
RPL=21

MINL=200
MAXR2=250
QUAL=1
```

## Pre-processing
Script will:<br>
1. Remove reads with both forward and reverse primers<br>
2. Remove reads with adapter contamination<br>
3. Filter for quality and minimum length (with UTRIM)<br>
4. Convert FASTQ to single line FASTA

```shell

$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$ARDERI/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $ARDERI/data/$RUN/$SSU \
 $ARDERI/metabarcoding_pipeline/primers/primers.db \
 $MINL $MAXR2 $QUAL; 
```

### SSU/58S/LSU removal 

#### Identify SSU, 5.8S  and LSU regions

This will create a large number of array jobs on the cluster

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
$ARDERI/data/$RUN/$SSU/fasta \
R1 \
$ARDERI/metabarcoding_pipeline/hmm/ssu_end.hmm \
$ARDERI/metabarcoding_pipeline/hmm/58s_start.hmm \
ssu 58ss 20

$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
 $ARDERI/data/$RUN/$SSU/fasta \
 R2 \
 $ARDERI/metabarcoding_pipeline/hmm/lsu_start.hmm \
 $ARDERI/metabarcoding_pipeline/hmm/58s_end.hmm \
 lsu 58se 20


```

#### Remove SSU, 5.8S  and LSU regions and merge output

If reverse read quality was poor and it was necessary to truncate reads to get more than a couple of reads past set LOWQUAL to TRUE

LOWQUAL keeps reads which lack 5.8S homology - this is necessary as trimming will in most instances have removed the homologous region

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$ARDERI/data/$RUN/$SSU/fasta/*R1" \
 $ARDERI/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58"

LOWQUAL=FALSE   
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$ARDERI/data/$RUN/$SSU/fasta/*R2" \
 $ARDERI/metabarcoding_pipeline/scripts/rm_58Se_LSU_v2.R \
 "*.\\.58" \
 "*.\\.lsu" \
 $LOWQUAL
```

Correct fasta if just R1 is to be used
```
for f in $ARDERI/data/$RUN/ITS/unfiltered/*r1*; do
F=$(echo $f|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}')
L=$(echo $f|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") 
awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $F.tmp && mv $F.tmp $F
done
```

#### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2

```shell
mkdir -p $ARDERI/data/$RUN/$SSU/filtered
find $ARDERI/data/$RUN/$SSU/fasta -type f -name *.r*.fa|xargs -I myfile mv myfile $ARDERI/data/$RUN/$SSU/filtered/.

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
FPL=23 
RPL=21

### Cluster 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $ARDERI $RUN $SSU 0 0
```

Work around for usearch bug 10.1
```shell
sed -ie 's/Zotu/OTU/' ITS.zotus.fa
```

### Assign taxonomy
```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $ARDERI $RUN $SSU 
```

### OTU evolutionary distance

Output a phylogentic tree in phylip format (both upper and lower triangles)
(usearch9 doesn't work)
```shell
usearch8.1 -calc_distmx ITS.otus.fa -distmxout ITS.phy -distmo fractdiff -format phylip_square
```

### Create OTU tables

Concatenates unfiltered reads, then assigns forward reads to OTUs. For any non-hits, attemps to assign reverse read (ITS2) to an OTU. 

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL true
```


### Combine biome and taxa

biom_maker.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biome and output a combined taxa and biome file to standard out.

```shell
$ARDERI/metabarcoding_pipeline/scripts/biom_maker.pl ITS.taxa ITS.otu_table.biom >ITS.taxa.biom
```


### [16S workflow](../master/16S%20%20workflow.md)

### [Statistical analysis](../master/statistical%20analysis.md)



## Old stuff

```shell


##### Concatenate
#cat $METAGENOMICS/data/$RUN/ITS/filtered/*.fa > $METAGENOMICS/data/$RUN/ITS.t.fa
##### Pad
#X=`cat ITS.t.fa|awk '{if ($1!~/>/) {print length($0)};}'|awk '$0>x{x=$0};END{print x}'`
#usearch8.1 -fastx_truncate ITS.t.fa -trunclen $X -padlen $X -fastaout ITS.fa
#rm ITS.t.fa
##### Dereplicate
#cat ITS.fa|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|$METAGENOMICS/scripts/get_uniq.pl > #ITS.sorted.fasta 
#rm ITS.fa
##### Cluster
#usearch8.1 -cluster_otus ITS.sorted.fasta -otus ITS.otus.fa -uparseout ITS.out.up -relabel OTU -minsize 2 
```
```shell
##### Concatenate unfiltered reads (Unfiltered fastq will need to be converted to fasta first)
for f in $METAGENOMICS/data/$RUN/ITS/unfiltered/*.r1.*
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	awk -v S="$S" -v FPL=$FPL -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,FPL);print $1}}' $f >>ITS1.unfiltered.fa
done

##### Make table (creates an OTU table of read counts per OTU per sample)
usearch9 -usearch_global ITS1.unfiltered.fa -db ITS.otus.fa -strand plus -id 0.97 -biomout ITS1.otu_table.biom -otutabout ITS1.otu_table.txt -output_no_hits -userout ITS1.hits.out -userfields query+target
```

```shell

for f in $METAGENOMICS/data/$RUN/ITS/unfiltered/*.r2.*
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	awk -v S="$S" -v RPL=$RPL -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,RPL);print $1}}' $f >t1
	grep "$S.*\*" ITS1.hits.out|awk -F";" '{print $2}'|awk -F" " '{print $1}'|$METAGENOMICS/scripts/seq_select_v2.pl t1 >> ITS2.unfiltered.fa
done
rm t1

usearch9 -usearch_global ITS2.unfiltered.fa -db ITS.otus.fa -strand both -id 0.97 -biomout ITS2.otu_table.biom -otutabout ITS2.otu_table.txt

```
