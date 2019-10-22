# folder containing project files
PROJECT_FOLDER=~/projects/Endophytes

# sequencer run folder 
RUN=leone_230919

# make the project folder
mkdir -p $PROJECT_FOLDER

# link to metabarcoding pipeline (specific for authors profile)
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

# folder to hold fatsq files
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq

# loop, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for s in "BAC FUN"; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
  mkdir $PROJECT_FOLDER/data/$RUN/$s/merged
done

# QC
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/slurm/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done

# Demultiplex bacterial and fungal amplicons
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/slurm/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*_R1_*" 0 \
 $P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.


# pre-process BAC files (min length 300, max diffs 5, quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/slurm/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 5 0.5 17 21

# Pre-process FUN files (min length 150, max diffs 5, quality 0.5) 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/slurm/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 200 1 23 21
 
   
