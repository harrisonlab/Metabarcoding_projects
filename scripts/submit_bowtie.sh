#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

R1=$1
R2=$2
REF=$3
DEST=$4
OUTNAME=$5
X=$7
I=$6
WORK_DIR=$TMPDIR

echo  "Running Bowtie 2S with the following in= REF IS '$REF' READ 1 '$R1' READ 2 ' $R2 ' DEST is '$DEST' OUTNAME IS '$OUTNAME' MIN is '$I' MAX is '$X'"

cd $WORK_DIR
#bowtie2 -p 4 --no-unal --un-conc $WORK_DIR -X $X -I $I -x $REF -1 $R1 -2 $R2 -S $WORK_DIR/$OUTNAME
#bowtie2 -p 8 --no-unal -X $X -I $I -x $REF -1 $R1 -2 $R2 -S $OUTNAME

###unmatched pairs
bowtie2 -p 8 --un-conc $OUTNAME -X $X -I $I -x $REF -1 $R1 -2 $R2 > /dev/null


cp * $DEST/.