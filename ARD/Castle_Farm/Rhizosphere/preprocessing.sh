# folder containing project files
PROJECT_FOLDER=~/projects/ARD

# sequencer run folder 
RUN=170217

# make the project folder
mkdir -p $PROJECT_FOLDER

# link to metabarcoding pipeline (specific for authors profile)
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

# folder to hold fatsq files
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq

# variable to hold folder names (BAC and FUN)
RIB="BAC FUN OO NEM"

# loop through the RIB variable, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for s in $RIB; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# QC
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/Kenya*; do
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

# Demultiplex nematode and oomycete amplicons
P1F=CGCGAATRGCTCATTACAACAGC
P1R=GGCGGTATCTGATCGCC
P2F=GAAGGTGAAGTCGTAACAAGG
P2R=AGCGTTCTTCATCGATGTGC

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*Nem*_R1_*" 0 \
 $P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/NEM/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/OO/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.


# pre-process BAC files (min length 300, max diffs 5, quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/BAC \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 5 0.5 17 21

# Pre-process FUN files (min length 200, MAX length 300, quality 1)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 200 1 23 21
# move FUN files to required location
for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
  FO=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
  L=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}') ;
  echo $L
  awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
done

# Pre-process OO files (min length 150, max diffs 10 (actual: (min len * max diffs)/100), quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OOpre \
 "$PROJECT_FOLDER/data/$RUN/OO/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/OO \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 150 10 0.5 21 20
 
# Pre-process NEM files (min length 100, max length, quality 1)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c NEMpre \
 "$PROJECT_FOLDER/data/$RUN/NEM/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/NEM \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/nematode.db \
 200 1 23 18
 # move FUN files to required location
for F in $PROJECT_FOLDER/data/$RUN/NEM/fasta/*_R1.fa; do 
  FO=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
  L=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}') ;
  echo $L
  awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/NEM/filtered/$FO;
done
