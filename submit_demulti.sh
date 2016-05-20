#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FORWARD=$1
REVERSE=$2
FPRIM1=$3
RPRIM1=$4
FPRIM2=$5
RPRIM2=$6
MISSMATCH=$7
S=$8

${S}/demulti.pl $FORWARD $REVERSE $FPRIM1 $RPRIM1 $FPRIM2 $RPRIM2 $MISSMATCH 


