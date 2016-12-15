#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

dir=$1

cat $OUTDIR/$dir/*.t1.fa > $dir/t1.fa
cat $OUTDIR/$dir/*.t2.fa > $dir/t2.fa


