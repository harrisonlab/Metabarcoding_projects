#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

OUTDIR=$1/data/$2
PREFIX=$3
UNFILTDIR=$OUTDIR/$PREFIX/unfiltered
SL=$4
SR=$5

for SCRIPT_DIR; do true; done

#### Concatenate unfiltered reads (Unfiltered fastq will need to be converted to fasta first )
for f in $UNFILTDIR/*.fastq
do
    S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
    R=$(echo $f|awk -F"." '{print $2}')
    if [ "$R" == "r1" ]; then
         $SCRIPT_DIR/fq2fa_v2.pl $f $OUTDIR/t1.fa $S $SL
    elif [ "$R" == "r2" ]; then
         $SCRIPT_DIR/fq2fa_v2.pl $f $OUTDIR/t2.fa $S 0 $SR
    else
        $SCRIPT_DIR/fq2fa_v2.pl $f $OUTDIR/t1.fa $S $SL $SR
    fi
done


#### Make table (Creates an OTU table of read counts per OTU per sample)

usearch9 -usearch_global t1.fa -db $PREFIX.otus.fa -strand plus -id 0.97 -biomout $PREFIX.otu_table.biom -otutabout $PREFIX.otu_table.txt -output_no_hits -userout $PREFIX.hits.out -userfields query+target

rm t1.fa
if [ "$R" == "r1" ] || [ "$R" == "r2" ];then
	for f in $UNFILTDIR/*.r2.*;do
		S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
		grep "$S.*\*" $PREFIX.hits.out|awk -F";" '{print $2}'|awk -F" " '{print $1}'|$SCRIPT_DIR/seq_select_v2.pl t2.fa >> t3.fa
	done
	rm t2.fa
	usearch9 -usearch_global t3.fa -db $PREFIX.otus.fa -strand both -id 0.97 -biomout ${PREFIX}2.otu_table.biom -otutabout ${PREFIX}2.otu_table.txt
	rm t3.fa
fi

rm $PREFIX.hits.out



