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
MAXR2=$7
QUAL=$8
SCRIPT_DIR=$9

mkdir -p $OUTDIR 
cd $OUTDIR 

usearch9 -search_oligodb $F -db $PRIMERS -strand both -userout ${F}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
echo "searched forward primers"

usearch9 -search_oligodb $R -db $PRIMERS -strand both -userout ${R}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
echo "searched reverse primer"


sed 's|^|/|;s|$|/,+3 d|' <(grep primer1 ${F}.t1.txt|awk -F"\t" '{print $1}') > ${OUTFILE}.temp.sed
echo "Created forward sed"
sed -f ${OUTFILE}.temp.sed $F > ${F}.t1.fastq
echo "searched forward"

sed 's|^|/|;s|$|/,+3 d|' <(grep primer2 ${R}.t1.txt|awk -F"\t" '{print $1}') > ${OUTFILE}.temp.sed
echo "search reverse sed"
sed -f ${OUTFILE}.temp.sed $R > ${R}.t1.fastq 
echo "searched reverse"

echo "deleting forward"
grep "adapter|primer2" -E ${F}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${F}.t1.fastq > ${F}.t2.fastq
echo "deleted forward, deleting reverse"
grep "adapter|primer1" -E ${R}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${R}.t1.fastq > ${R}.t2.fastq
echo "deleted reverse"


usearch9 -fastq_filter ${F}.t2.fastq -fastq_minlen $MINL -fastq_maxee $QUAL -fastaout ${OUTFILE}_t1.fa
usearch9 -fastq_filter ${R}.t2.fastq -fastq_minlen $MINL -fastq_trunclen $MAXR2 -fastq_maxee $QUAL -fastaout ${OUTFILE}_t2.fa

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
