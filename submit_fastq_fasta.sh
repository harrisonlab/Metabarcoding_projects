#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

INFILE=$(sed -n -e "$SGE_TASK_ID p" $1)
dir=$2
SL=$3
SR=$4
FILE=$(echo $INFILE|awk -F"/" '{print $NF}')
PREFIX=$(echo $FILE|awk -F"." '{print $1}')

for SCRIPT_DIR; do true; done

if [[ "$FILE" =~ "r1" ]]; then
	$SCRIPT_DIR/fq2fa_v2.pl $INFILE $dir/$PREFIX.t1.fa $PREFIX $SL 0
elif [[ "$FILE" =~ "r2" ]]; then 
	$SCRIPT_DIR/fq2fa_v2.pl $INFILE $dir/$PREFIX.t2.fa $PREFIX 0 $SR
else
	$SCRIPT_DIR/fq2fa_v2.pl $INFILE $dir/$PREFIX.t1.fa $PREFIX $SL $SR
fi
