# Common workflow

Set ARDERI to project folder.

RUN should be set to location where files are to be stored.

```shell
#run for each analysis
mkdir -p $ARDERI/data/$RUN/fastq
mkdir $ARDERI/data/$RUN/quality
mkdir $ARDERI/data/$RUN/ambiguous
mkdir -p $ARDERI/data/$RUN/16S/fastq
mkdir $ARDERI/data/$RUN/16S/filtered
mkdir $ARDERI/data/$RUN/16S/unfiltered
mkdir -p $ARDERI/data/$RUN/ITS/fastq
mkdir $ARDERI/data/$RUN/ITS/filtered
mkdir $ARDERI/data/$RUN/ITS/unfiltered
```

## Decompress files
```shell
for FILE in $ARDERI/data/$RUN/fastq/*.gz; do 
	$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c unzip $FILE
done
```

## QC
Qualtiy checking with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```shell
for FILE in $ARDERI/data/$RUN/fastq/*; do 
	$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $ARDERI/data/$RUN/quality
done
```

## Demultiplexing

We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Ambiguous reads are written to two (f & r) seperate files.

Running something like the below should give a good indication of what index_1 and index_2 should be - this is useful if you don't knwo what the primer sequences are and to get a feel of how many mismatches (if necesary) to use. 
```shell
sed -n '2~4p' $(ls|head -n1)|grep -x "[ATCG]\+"|cut -c-16|sort|uniq| \
tee zzexpressions.txt|xargs -I%  grep -c "^%" $(ls|head -n1) >zzcounts.txt
```

Any sequence which has too many mismatches, or none mathching primers is removed to a file x.ambigous.fq

demultiplex can accept any number of primer pairs (though for this project only 2 primer pairs are multiplexed)

<table>
Primers:
<tr><td><td>Forward<td>Reverse</tr>
<tr><td>16S<td>CCTACGGGNGGCWGCAG<td>GACTACHVGGGTATCTAATCC</tr>
<tr><td>ITS<td>CTTGGTCATTTAGAGGAAGTAA<td>ATATGCTTAAGTTCAGCGGG</tr>
<tr><td>OO<td>GAAGGTGAAGTCGTAACAAGG<td>AGCGTTCTTCATCGATGTGC</tr>
<tr><td>Nem<td>CGCGAATRGCTCATTACAACAGC<td>GGCGGTATCTGATCGCC</tr>
</table>


```shell
P1F=CCTACGGGNGGCWGCAG # bacteria
P1R=GACTACHVGGGTATCTAATCC 
P2F=CTTGGTCATTTAGAGGAAGTAA # fungi
P2R=ATATGCTTAAGTTCAGCGGG

$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
	"$ARDERI/data/$RUN/fastq/*_R1_*" 0 \
	$P1F $P1R $P2F $P2R


mv $ARDERI/data/$RUN/fastq/*ps1* $ARDERI/data/$RUN/16S/fastq/.
mv $ARDERI/data/$RUN/fastq/*ps2* $ARDERI/data/$RUN/ITS/fastq/.
mv $ARDERI/data/$RUN/fastq/*ambig* $ARDERI/data/$RUN/ambiguous/.
```

###[16S workflow](../master/16S%20%20workflow.md)
###[ITS workflow](../master//ITS%20workflow.md)
