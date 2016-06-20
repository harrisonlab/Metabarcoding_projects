#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
ADAPTERS=$5
MINL=$6
QUAL=$7
LABEL=$8

mkdir -p $OUTDIR 
cd $OUTDIR 

echo usearch8.1 -fastq_mergepairs $F -reverse $R -fastqout ${OUTFILE}.t1 
echo usearch8.1 -search_oligodb ${OUTFILE}.t1 -db $ADAPTERS -strand both -userout ${OUTFILE}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
echo cat ${OUTFILE}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I � sed -i -ne:t -e"/*\@�.*/D" -e'$!N;//D;/'"\@�/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' ${OUTFILE}.t1
echo usearch8.1 -fastq_filter ${OUTFILE}.t1 -fastq_minlen $MINL -fastq_maxee_rate $QUAL -relabel $LABEL -fastaout $OUTFILE

echo mkdir -p $OUTDIR/../unfiltered 
echo mv ${OUTFILE}.t1 $OUTDIR/../unfiltered/${OUTFILE}.unfiltered.fastq

echo rm ${OUTFILE}.t1.txt