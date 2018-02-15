# folder containing project files
PROJECT_FOLDER=~/projects/metabarcoding/ARD

# sequencer run folders 
RUNS="180105 180109 180112 180119"

# common folders
for RUN in $RUNS; do
  #mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/quality
  mkdir $PROJECT_FOLDER/data/$RUN/ambiguous
done

# variable to hold folder names (BAC and FUN)
RIB="BAC FUN OO NEM"

# loop through the RIB variable, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for RUN in $RUNS; do
  for s in $RIB; do
    mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
    mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
    mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
    mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
  done
done

# QC
for RUN in $RUNS; do
  for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do
    $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
  done
done

# BAC and FUN are multiplexed. Can seperate by the primer sequences (p1 for BAC, p2 for FUN)
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

# demultiplex with 0 difference in primer seqeunce
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
   "$PROJECT_FOLDER/data/$RUN/fastq/*R1*.gz" 0 \
   $P1F $P1R $P2F $P2R
done

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/16S/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/ITS/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.

# pre-process 16S files (min length 300, max diffs 5, quality 0.5)
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
   "$PROJECT_FOLDER/data/$RUN/16S/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/16S \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
   300 5 0.5 0 0
done

# Pre-process FUN files (min length 200, MAX length 300, quality 1)
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
   "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/ITS \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
   200 1 0 0
done
