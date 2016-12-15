#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

OUTDIR=$1/data/$2
PREFIX=$3

for SCRIPT_DIR; do true; done

cd $OUTDIR

#### Assign Taxonomy
usearch9 -utax ${PREFIX}.otus.fa -db $SCRIPT_DIR/../taxonomies/utax/${PREFIX}_ref.udb -strand both -utaxout ${PREFIX}.reads.utax -rdpout ${PREFIX}.rdp -alnout ${PREFIX}.aln.txt
cat ${PREFIX}.rdp|$SCRIPT_DIR/mod_taxa.pl > ${PREFIX}.taxa 
