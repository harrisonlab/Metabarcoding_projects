#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FILTDIR=$1
OUTDIR=$2
PREFIX=$3
SL=$4
SR=$5

for SCRIPT_DIR; do true; done

cd $OUTDIR

#### Concatenate files
cat ${FILTDIR}/*.fa > ${PREFIX}.temp.fa

#### Remove multiplex primers and pad reads to same length
X=`awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}' ${PREFIX}.temp.fa|awk '$0>x{x=$0};END{print x}'`
usearch8.1 -fastx_truncate ${PREFIX}.temp.fa -stripleft $SL -stripright $SR -trunclen $X -padlen $X -fastaout ${PREFIX}.fa
rm ${PREFIX}.temp.fa

#### Dereplication
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${PREFIX}.fa|$SCRIPT_DIR/get_uniq.pl > ${PREFIX}.sorted.fasta 
rm ${PREFIX}.fa

#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))
usearch8.1 -cluster_otus ${PREFIX}.sorted.fasta -otus ${PREFIX}.otus.fa -uparseout ${PREFIX}.out.up -relabel OTU -minsize 2
rm ${PREFIX}.sorted.fasta

#### Assign Taxonomy
#usearch8.1 -utax ${PREFIX}.otus.fa -db $SCRIPT_DIR/../taxonomies/utax/${PREFIX}_ref.udb -strand both -utaxout ${PREFIX}.reads.utax -rdpout ${PREFIX}.rdp -alnout ${PREFIX}.aln.txt
#cat ${PREFIX}.rdp|$SCRIPT_DIR/mod_taxa.pl > ${PREFIX}.taxa 