# make project folders
PROJECT_FOLDER=~/projects/Ash
mkdir -p $PROJECT_FOLDER
ln -s $MBPL $PROJECT_FOLDER/metabarcoding_pipeline 

RUN=.
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous
mkdir $PROJECT_FOLDER/data/$RUN/cluster

for s in BAC FUN; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# quality check
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do 
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/slurm/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done

# Demultiplex bacterial and fungal amplicons
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*_R1_*" 0 \
 $P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.

