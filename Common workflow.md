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

We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Ambiguous reads are written to two (f & r) seperate files.

Running something like the below should give a good indication of what index_1 and index_2 should be - this is useful if you don't knwo what the primer sequences are and to get a feel of how many mismatches (if necesary) to use. 
```shell
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-16|sort|uniq > zzexpressions.txt
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-16|sort|uniq|xargs -I r grep -c ^r $(ls|head -n1) >zzcounts.txt
```

I typically use the first 8 nucleotides of the primer and allow 0 mismatches (the final parameter). Any sequence which has too many mismatches, or none mathching primers is removed to a file x.ambigous.fq

demultiplex can accept any number of primer pairs (though for this project only 2 primer pairs are multiplexed)

```shell
P1F=CCTACGGG # bacteria
P1R=GACTACHV 
P2F=CTTGGTCA # fungi
P2R=ATATGCTT

METAGENOMICS/scripts/ARDERI.sh -c demultiplex /
	$METAGENOMICS/data/$RUN/fastq/*_R1_* 0/
	$P1F $P1R $P2F $P2R
```
	Type	F	R
	16S	CCTACGGG	GACTACHV
	ITS	CTTGGTCA	ATATGCTT
	ITS4-6	GAAGGTGA	TCCTCCGC
	ITS7-6	GAAGGTGA	AGCGTTCT
	Nem	CGCGAATR	GGCGGTAT

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
