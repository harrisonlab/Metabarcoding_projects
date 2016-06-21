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
SCRIPT_DIR=$9

LABEL=${OUTFILE}.

mkdir -p $OUTDIR 
cd $OUTDIR 

usearch8.1 -search_oligodb $F -db $PRIMERS -strand both -userout ${F}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
usearch8.1 -search_oligodb $R -db $PRIMERS -strand both -userout ${R}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

sed 's|^|/|;s|$|/,+3 d|' <(grep primer1 ${F}.t1.txt|awk -F"\t" '{print $1}') > ${OUTFILE}.temp.sed
sed -f ${OUTFILE}.temp.sed $F > ${F}.t1.fastq
sed 's|^|/|;s|$|/,+3 d|' <(grep primer2 ${R}.t1.txt|awk -F"\t" '{print $1}') > ${OUTFILE}.temp.sed
sed -f ${OUTFILE}.temp.sed $R > ${R}.t1.fastq 

grep adapter ${F}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${F}.t1.fastq > ${F}.t2.fastq
grep adapter ${R}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${R}.t1.fastq > ${R}.t2.fastq


usearch8.1 -fastq_filter ${F}.t2.fastq -fastq_minlen $MINL -fastq_maxee_rate $QUALF -relabel $LABEL -fastaout ${OUTFILE}_t1.fa
usearch8.1 -fastq_filter ${R}.t2.fastq -fastq_minlen $MINL -fastq_maxee_rate $QUALR -relabel $LABEL -fastaout ${OUTFILE}_t2.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}_t1.fa > ${OUTFILE}_R1.fa
sed -i -e '1d' ${OUTFILE}_R1.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}_t2.fa > ${OUTFILE}_R2.fa
sed -i -e '1d' ${OUTFILE}_R2.fa

mkdir -p $OUTDIR/../unfiltered 
mv ${F}.t2.fastq $OUTDIR/../unfiltered/${OUTFILE}.r1.unfiltered.fastq
mv ${R}.t2.fastq $OUTDIR/../unfiltered/${OUTFILE}.r2.unfiltered.fastq


rm ${F}.t1.txt ${R}.t1.txt ${F}.t1.fastq ${R}.t1.fastq ${OUTFILE}_t1.fa ${OUTFILE}_t2.fa ${OUTFILE}.temp.sed


#xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' ${F}.t.fastq
#xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' ${R}.t.fastq
