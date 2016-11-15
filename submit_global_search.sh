#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

INFILE=$1
OUTDIR=$2
PREFIX=$3
EP=$4


if [ -z $EP ]; then
	usearch9 -usearch_global $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand plus -id 0.97 -biomout $OUTDIR/$PREFIX.otu_table.biom -otutabout $OUTDIR/$PREFIX.otu_table.txt -output_no_hits -userout $OUTDIR/$PREFIX.hits.out -userfields query+target
else
	usearch9 -usearch_global $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand both -id 0.97 -biomout $OUTDIR/$PREFIX$EP.otu_table.biom -otutabout $OUTDIR/$PREFIX$EP.otu_table.txt -output_no_hits -userout $OUTDIR/$PREFIX$EP.hits.out -userfields query+target
fi