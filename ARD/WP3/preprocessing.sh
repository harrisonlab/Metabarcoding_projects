YEAR1="151019 160408 160418 160506 160111 160122 160201" 
YEAR2="170110 170224 170323 170424"
YEAR3="180105 180109 180112 180119"

RUNS="$YEAR1 $YEAR2 $YEAR3"

for RUN in $RUNS; do
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
for RUN in $RUNS; do
  for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do 
    $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
  done
done

# demultiplex
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
   "$PROJECT_FOLDER/data/$RUN/fastq/*16s*_R1_*" 0 \
   $P1F $P1R $P2F $P2R
done

for RUN in $RUNS; do
  mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
  mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
  mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.
done

# BACTERIA
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
   "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/BAC \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
   300 5 0.5 17 21 
done 

 # FUNGI
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
   "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/FUN \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
   200 1 23 21
done

for RUN in $RUNS; do
  for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
   FO=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1".r1.fa"}'); 
   L=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1}') ;
   echo $L
   awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
  done
done

# Oomycetes
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OOpre \
   "$PROJECT_FOLDER/data/$RUN/OO/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/OO \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
   150 5 0.5 21 20
done    

# Nematode
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c NEMpre \
   "$PROJECT_FOLDER/data/$RUN/NEM/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/NEM \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/nematode.db \
   150 10 0.5 23 18
done


# make analysis folders
mkdir $PROJECT_FOLDER/analysis/WP3
for s in BAC FUN OO NEM; do
  mkdir -p $PROJECT_FOLDER/analysis/WP3/$s/filtered
  mkdir $PROJECT_FOLDER/analysis/WP3/$s/unfiltered
done

# link processed data to anaylsis folders
for RUN in $RUNS; do
  for s in BAC FUN; do
    ln -s $PROJECT_FOLDER/data/$RUN/$s/unfiltered/* $PROJECT_FOLDER/analysis/WP3/$s/unfiltered/.
    ln -s $PROJECT_FOLDER/data/$RUN/$s/filtered/* $PROJECT_FOLDER/analysis/WP3/$s/filtered/.
    ln -s $PROJECT_FOLDER/data/$RUN/ambiguous/* $PROJECT_FOLDER/analysis/WP3/ambiguous/.
  done
done
