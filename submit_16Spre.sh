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
SCRIPT_DIR=$8

LABEL=${OUTFILE}.

mkdir -p $OUTDIR 
cd $OUTDIR 

usearch8.1 -fastq_mergepairs $F -reverse $R -fastqout ${OUTFILE}.t1  -fastq_maxdiffpct 15 -fastq_maxdiffs $(($MINL*15/100)) 
usearch8.1 -search_oligodb ${OUTFILE}.t1 -db $ADAPTERS -strand both -userout ${OUTFILE}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 


cat ${OUTFILE}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${OUTFILE}.t1 > ${OUTFILE}.t2

#xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' ${OUTFILE}.t1

usearch8.1 -fastq_filter ${OUTFILE}.t2 -fastq_minlen $MINL -fastq_maxee $QUAL -relabel $LABEL -fastaout ${OUTFILE}.t3.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}.t3.fa > ${OUTFILE}.filtered.fa
sed -i -e '1d' ${OUTFILE}.filtered.fa


mkdir -p $OUTDIR/../unfiltered 
mv ${OUTFILE}.t2 $OUTDIR/../unfiltered/${OUTFILE}.unfiltered.fastq

rm ${OUTFILE}.t1.txt ${OUTFILE}.t1 ${OUTFILE}.t3.fa
