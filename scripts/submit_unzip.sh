#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FILE=$1

for SCRIPT_DIR; do true; done

pigz -d $FILE

