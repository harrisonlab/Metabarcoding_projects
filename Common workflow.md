# Common workflow

```shell
#run for each analysis
mkdir -p $METAGENOMICS/analysis/$RUN/16S
mkdir $METAGENOMICS/analysis/$RUN/ITS	
mkdir -p $METAGENOMICS/data/$RUN/fastq
mkdir -p $METAGENOMICS/data/$RUN/16S/fastq
mkdir $METAGENOMICS/data/$RUN/16S/filtered
mkdir $METAGENOMICS/data/$RUN/16S/unfiltered
mkdir -p $METAGENOMICS/data/$RUN/ITS/fastq
mkdir $METAGENOMICS/data/$RUN/ITS/filtered
mkdir $METAGENOMICS/data/$RUN/ITS/unfilterd

```	
The $METAGENOMICS directory should be set to something appropriate (e.g. /home/bob/metagenomics) and $RUN to the name of the NGS run. The realtive path is used in the scripts below - depending on your config you may have to specify full paths.	


## QC
Qualtiy checking was performed with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

From same folder containing fastq files ran:

	fastqc *

## Demultiplexing

### OLD (fast) method
We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Ambiguous reads are written to two (f & r) seperate files.

Running something like the below should give a good indication of what index_1 and index_2 should be - this is useful if you don't knwo what the primer sequences are and to get a feel of how many mismatches to use. 
```shell
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-16|sort|uniq > zzexpressions.txt
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-16|sort|uniq|xargs -I r grep -c ^r $(ls|head -n1) >zzcounts.txt
```

I typically use the first 8 nucleotides of the primer and allow 2 mismatches (the final parameter). Anything sequence which has too many mismatches, or none mathching primers is removed to a file x.ambigous.fq
```shell
for f in $METAGENOMICS/data/$RUN/fastq/*_R1_*
do     
	R1=$f     
	R2=$(echo $R1|sed 's/_R1_/_R2_/')    
	S=$(echo $f|awk -F"_" '{print $2}')     
	echo $f    
	$METAGENOMICS/scripts/demulti.sh $R1 $R2 "CCTACGGG" "GACTACHV" "CTTGGTCA" "ATATGCTT" 2  
done   
```
	F	R		
ITS4-6	TCCTCCGC	GAAGGTGA
Nem	CGCGAATR	GGCGGTCC	
ITS7-6	GAAGGTGA	AGCGTTCT


```shell
mkdir -p $METAGENOMICS/data/$RUN/16S/fastq
mkdir -p $METAGENOMICS/data/$RUN/ITS/fastq
mkdir -p $METAGENOMICS/data/$RUN/ambiguous

mv $METAGENOMICS/data/$RUN/fastq/*bacterial* $METAGENOMICS/data/$RUN/16S/fastq/.
mv $METAGENOMICS/data/$RUN/fastq/*fungal* $METAGENOMICS/data/$RUN/ITS/fastq/.
mv $METAGENOMICS/data/$RUN/fastq/*ambig* $METAGENOMICS/data/$RUN/ambiguous/.
```

###[16S workflow](../master/16S%20%20workflow.md)
###[ITS workflow](../master//ITS%20workflow.md)
