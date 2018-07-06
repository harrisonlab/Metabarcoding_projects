# folder containing project files
PROJECT_FOLDER=~/projects/Oak_decline

# sequencer run folder 
RUN=130618

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
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done

# no multiplexing but demulti_v3.pl has some useful features for discovering barcodes per sample and filtering on barcode + primer errors
# BAC discover sample barcode + filter with max 1 difference barcode + primer (forward or reverse)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*R1*.gz" 11 \
 CCTAYGGGNGGCWGCAG GGACTACNNGGGTATCTAATCC 
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.

# FUN discover sample barcode + filter with max 1 difference barcode + primer (forward or reverse)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
 "$PROJECT_FOLDER/data/$RUN/fastq/*R1*.gz" 11 \
 GGAAGTAAAAGTCGTAACAAGG GCTGCGTTCTTCATCGATGC 
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.


# pre-process BAC files (min length 150, max diffs 5, quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/BAC \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/all_primers_dapters.db \
 150 5 0.5 23 28

# Pre-process FUN files (min length 150, max diffs 5, quality 0.5) 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/all_primers_dapters.db \
 150 5 0.5 28 26
 
# remove unfiltered reads - reads without barcode and primer have been supplied
rm  $PROJECT_FOLDER/data/$RUN/BAC/unfiltered/*
rm  $PROJECT_FOLDER/data/$RUN/FUN/unfiltered/*

# merge pre-cut reads
for s in "BAC FUN"; do   
  for FORWARD in $PROJECT_FOLDER/data/$RUN/$s/unmerged/*_1.fq; do
    REVERSE=$(sed 's/_1.fq/_2.fq/' <<< $FORWARD) 
    OUTFILE=$(sed 's/_1.fq//' <<< $FORWARD|sed 's/.*\///')
    qsub $PROJECT_FOLDER/metabarcoding_pipeline/scripts/submit_merge_only.sh \
    $FORWARD $REVERSE $OUTFILE $PROJECT_FOLDER/data/$RUN/$s/merged 150 5 
  done  
done
   
