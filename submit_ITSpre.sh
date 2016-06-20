#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
PRIMERS=$5
MINL=$6
QUALF=$7
QUALR=$8
LABEL=$9

mkdir -p $OUTDIR 
cd $OUTDIR 

usearch8.1 -search_oligodb $F -db $PRIMERS -strand both -userout ${F}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
usearch8.1 -search_oligodb $R -db $PRIMERS -strand both -userout ${R}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 


sed 's|^|/|;s|$|/,+3 d|' <(grep primer1 ${F}.t1.txt|awk -F"\t" '{print $1}') > temp.sed
sed -f temp.sed $F > ${F}.t.fastq
sed 's|^|/|;s|$|/,+3 d|' <(grep primer2 ${R}.t1.txt|awk -F"\t" '{print $1}') > temp.sed
sed -f temp.sed $F > ${R}.t.fastq 


grep adapter ${F}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' ${F}.t.fastq
grep adapter ${R}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' ${R}.t.fastq


usearch8.1 -fastq_filter ${F}.cleaned.fastq -fastq_minlen $MINL -fastq_maxee_rate $QUALF -relabel $LABEL -fastaout $OUTFILE
usearch8.1 -fastq_filter ${R}.cleaned.fastq -fastq_minlen $MINL -fastq_maxee_rate $QUALR -relabel $LABEL -fastaout $OUTFILE


rm ${F}.t1.txt ${R}.t1.txt ${F}.t.fastq ${R}.t.fastq
