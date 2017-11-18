# folder containing project files
PROJECT_FOLDER=~/projects/ARD

# sequence data location(s)
RUN=170110
RUN=170224
RUN=170323
RUN=170424 #nem/oo

# make the project folder
mkdir -p $PROJECT_FOLDER

# link to metabarcoding pipeline (specific for authors profile)
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

# folder to hold fastq files
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq

# variable to hold folder names (BAC and FUN)
RIB="BAC FUN"

# loop through the RIB variable, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for s in $RIB; do
mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# QC
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*.gz; do
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done

# BAC and FUN are multiplexed. Can seperate by the primer sequences (p1 for BAC, p2 for FUN)
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

# demultiplex with 0 difference in primer seqeunce
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
"$PROJECT_FOLDER/data/$RUN/fastq/*R1*.gz" 0 \
$P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.

# bacteria
SSU=BAC;FPL=17;RPL=21
MINL=300;MINOVER=5;QUAL=0.5

## pre-process BAC files (min length 300, max diffs 5, quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
$PROJECT_FOLDER/data/$RUN/$SSU \
$PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL

# fungi
SSU=FUN;FPL=23;RPL=21
MINL=200;MAXR2=300;QUAL=1

## Pre-process FUN files (min length 200, MAX R2 length 300, quality 1)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 $MINL $MAXR2 $QUAL; 

for F in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*_R1.fa; do 
 FO=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}') ;
 echo $F
 echo $FO
 echo $L
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$FO;
done


##### BELOW NO LONGER IMPLEMENTED
## identify none ITS (FUN) regions (R1 only)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
$PROJECT_FOLDER/data/$RUN/$SSU/fasta \
R1 \
$PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
$PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
ssu 58ss 20

### remove none ITS regions from sequence
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fasta/*R1" \
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58" \
 $FPL $RPL
 
### Move merged fasta to filtered folder
for D in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*1; do 
 F=$(echo $D|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $D|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") ;
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $D/$F > $F.tmp && mv $F.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$F;
# mv $PROJECT_FOLDER/data/$RUN/$SSU/fasta/${L}_R1/$F $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$L; 
done

# move files
RUN="170110 170224 170323"
SSU="BAC FUN"
LOC="Goatham Heineken"

for R in $RUN; do
 for S in $SSU; do
  for L in $LOC; do
   T="${L:0:1}"
    ln -s ~/projects/ARD/data/$R/$S/filtered/${T}1* ~/projects/ARD/analysis/$L/Year2/$S/filtered/.
    ln -s ~/projects/ARD/data/$R/$S/unfiltered/${T}1* ~/projects/ARD/analysis/$L/Year2/$S/unfiltered/.
  done
 done
done




