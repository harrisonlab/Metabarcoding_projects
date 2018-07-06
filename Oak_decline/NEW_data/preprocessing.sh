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
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 150 5 0.5 23 28

# Pre-process FUN files (min length 150, max diffs 5, quality 0.5) 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 150 5 0.5 23 28
 
# move FUN files to required location
for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
  FO=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1".r1.fa"}'); 
  L=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1}') ;
  echo $L
  awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
done

# Pre-process OO files (min length 150, max diffs 10 (actual: (min len * max diffs)/100), quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OOpre \
 "$PROJECT_FOLDER/data/$RUN/OO/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/OO \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 150 10 0.5 21 20
 
# Pre-process NEM files (min length 150, max diffs 10 (actual: (min len * max diffs)/100), quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c NEMpre \
 "$PROJECT_FOLDER/data/$RUN/NEM/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/NEM \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/nematode.db \
 150 10 0.5 23 18
# move NEM files to required location
for F in $PROJECT_FOLDER/data/$RUN/NEM/fasta/*_R1.fa; do 
  FO=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1".r1.fa"}'); 
  L=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1}') ;
  echo $L
  awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/NEM/filtered/$FO;
done
