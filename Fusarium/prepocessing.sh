# make project folders
PROJECT_FOLDER=~/projects/fusarium
mkdir -p $PROJECT_FOLDER
ln -s $MBPL $PROJECT_FOLDER/metabarcoding_pipeline 

RUN=set_1
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous
mkdir $PROJECT_FOLDER/data/$RUN/cluster

for s in BAC FUN OG1 OG4; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# quality check
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do 
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done

# Demultiplex primer pairs
# 16S
#P0=16S
P0F=CCTACGGG[ATGC]GGC[AT]GCAG
P0R=GACTAC[ACT][ACG]GGGTATCTAATCC
#ITS
#P1=ITS
P1F=GTGAATCATCGAATCTTTGAACGC
P1R=CCGCTTATTGATATGCTTAA[AG]TTCAG
# TEF
#P2=TEF
P2F=GGTCACTTGATCTACCAGTGCG
P2R=CCCA[AG]GCGTACTTGAAG[AG]AAC
# SIX5
#P3=SIX5
P3F=AATCATCCTCACGATTATGTGTG
P3R=CTGATGGCAAAGGTCATAGAATGTT
# T4 orthogroup 13890
#P5=OG13890
P5F=GCTGTCTTATCACTTATCAGCCTTG
P5R=CGGTCTGATTTGGTGTCCAGTCG
# T6 orthogroup OG4952
#P6=OG4952
P6F=CCACACTTGACATGAGGAT[AG]GTC
P6R=GCTCACGGTCAGATAACTTTGC


#P1F=CCTACGGGNGGCWGCAG
#P1R=GACTACHVGGGTATCTAATCC
#P2F=CTTGGTCATTTAGAGGAAGTAA
#P2R=ATATGCTTAAGTTCAGCGGG

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*_R1_*" 0 \
 $P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.

# pre-process BAC files (min length 300, max diffs 5, quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/BAC \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 5 0.5 17 21 400

# Pre-process FUN files (min length 150, quality 1, truncate final 5 bases) 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 150 1 22 21 5
 
# move fasta files 
for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
 FO=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1}') ;
 echo $L
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
done   
