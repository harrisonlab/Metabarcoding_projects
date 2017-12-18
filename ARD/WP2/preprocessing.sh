PROJECT_FOLDER=~/projects/ARD
RUNS="171129 171215"

# make project folders
for RUN in RUNS; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/quality
  mkdir $PROJECT_FOLDER/data/$RUN/ambiguous
  for s in BAC FUN OO NEM; do
    mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
    mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
    mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
    mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
  done
done

# quality check
for RUN in RUNS; do
  for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do 
    $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
  done
done

# demultiplex
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

for RUN in RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
   "$PROJECT_FOLDER/data/$RUN/fastq/*16s*_R1_*" 0 \
   $P1F $P1R $P2F $P2R
done

for RUN in RUNS; do
  mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
  mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
  mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.
done

# BACTERIA
for RUN in RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
   "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/BAC \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
   300 5 0.5 17 21 
done 

 # FUNGI
 for RUN in RUNS; do
   $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
    "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
    $PROJECT_FOLDER/data/$RUN/FUN \
    $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
    200 1 23 21
done
