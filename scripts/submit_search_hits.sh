#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

INFILE=$(sed -n -e "$SGE_TASK_ID p" $1)
OUTDIR=$2
HITS=$3
S=$(echo $INFILE|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')

for SCRIPT_DIR; do true; done

grep "$S.*\*" $HITS|awk -F";" '{print $2}'|awk -F" " '{print $1}'|$SCRIPT_DIR/seq_select_v2.pl $OUTDIR/t2.fa >> $OUTDIR/t3.fa

