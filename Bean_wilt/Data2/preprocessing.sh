# folder containing project files
PROJECT_FOLDER=~/projects/Kenya

# sequencer run folder (in this instance all comparable data was in the same run)
RUN=Data2

# folder to hold fatsq files
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq

# variable to hold folder names (BAC and FUN)
RIB="FUN OO"

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

# demultiplex oo and fun
P1F=GAAGGTGAAGTCGTAACAAGG
P1R=AGCGTTCTTCATCGATGTGC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

# demultiplex with 0 difference in primer seqeunce
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
"$PROJECT_FOLDER/data/$RUN/fastq/*" 1 \
$P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/OO/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/. 

 # FUNGI
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
  "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
  $PROJECT_FOLDER/data/$RUN/FUN \
  $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
  200 1 23 21

 for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
  FO=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1".r1.fa"}'); 
  L=$(awk -F"/" '{print $NF}' <<< $F|awk -F"_" '{print $1}') ;
  echo $L
  awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
 done

# Oomycetes
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OOpre \
 "$PROJECT_FOLDER/data/$RUN/OO/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/OO \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 150 5 0.5 21 20
  
  
# make analysis folders
mkdir $PROJECT_FOLDER/analysis/$RUN
for s in FUN OO; do
  mkdir -p $PROJECT_FOLDER/analysis/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/analysis/$RUN/$s/unfiltered
done

# link processed data to anaylsis folders
for s in FUN OO; do
  ln -s $PROJECT_FOLDER/data/$RUN/$s/unfiltered/* $PROJECT_FOLDER/analysis/$RUN/$s/unfiltered/.
  ln -s $PROJECT_FOLDER/data/$RUN/$s/filtered/* $PROJECT_FOLDER/analysis/$RUN/$s/filtered/.
done
