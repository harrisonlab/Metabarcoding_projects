# make project folders
PROJECT_FOLDER=~/projects/fusarium
mkdir -p $PROJECT_FOLDER
ln -s $MBPL $PROJECT_FOLDER/metabarcoding_pipeline 

RUN=set_1
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous
mkdir $PROJECT_FOLDER/data/$RUN/cluster

for s in TEF OG1 OG4 SIX; do
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
P0F=CCTACGGGNGGCWGCAG
P0R=GACTACHVGGGTATCTAATCC
#ITS
#P1=ITS
P1F=GTGAATCATCGAATCTTTGAACGC
P1R=CCGCTTATTGATATGCTTAARTTCAG
# TEF
#P2=TEF
P2F=GGTCACTTGATCTACCAGTGCG
P2R=CCCARGCGTACTTGAAGRAAC
# SIX
#P3=SIX
P3F=AATCATCCTCACGATTATGTGTG
P3R=CTGATGGCAAAGGTCATAGAATGTT
# T4 orthogroup 13890
#P5=OG13890
P5F=GCTGTCTTATCACTTATCAGCCTTG
P5R=CGGTCTGATTTGGTGTCCAGTCG
# T6 orthogroup OG4952
#P6=OG4952
P6F=CCACACTTGACATGAGGATRGTC
P6R=GCTCACGGTCAGATAACTTTGC


#P1F=CCTACGGGNGGCWGCAG
#P1R=GACTACHVGGGTATCTAATCC
#P2F=CTTGGTCATTTAGAGGAAGTAA
#P2R=ATATGCTTAAGTTCAGCGGG
# s1 = P2, pS5?, ps6,foc= p2,p3?,p5,p6

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*_R1_*" 1 \
 $P2F $P2R $P5F $P5R $P6F $P6R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/TEF/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/SIX/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps3* $PROJECT_FOLDER/data/$RUN/OG1/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps4* $PROJECT_FOLDER/data/$RUN/OG4/fastq/.

mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.


# pre-process OG1 files (min length 300, max diffs 10 (30), quality 0.5,min merged length 356)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/OG1/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/OG1 \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 10 0.5 25 23 356

# pre-process OG4 files (min length 250, max diffs 12 (30), quality 0.5,min merged length 356)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/OG4/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/OG4 \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 250 12 0.5 23 22 277

# pre-process TEF files (min length 250, max diffs 12 (30), quality 0.5,min merged length 356)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/OG4/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/OG4 \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 250 12 0.5 23 22 320
 

# move fasta files 
for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
 FO=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1}') ;
 echo $L
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
done   
