PROJECT_FOLDER=~/projects/Cotton/
mkdir -p $PROJECT_FOLDER
ln -s $MBPL $PROJECT_FOLDER/metabarcoding_pipeline 

RUN=.

for s in BAC FUN; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

#bacteria
SSU=BAC
FPL=0
RPL=0

MINL=100
MINOVER=5
QUAL=0.5

# note PIPELINE.sh is set to use a different R1/R2 naming convention to that used with these files
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*_1.fq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 $MINL $MINOVER $QUAL $FPL $RPL 

#Fungi
SSU=FUN
FPL=0 
RPL=0

MINL=100
MINOVER=5
QUAL=0.5

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*_1.fq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 $MINL $MINOVER $QUAL $FPL $RPL 


RUN=RHIZOSPHERE
for s in BAC FUN; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mv $PROJECT_FOLDER/data/$s/filtered/*R* $PROJECT_FOLDER/data/$RUN/$s/filtered/. 
  mv $PROJECT_FOLDER/data/$s/unfiltered/*R* $PROJECT_FOLDER/data/$RUN/$s/unfiltered/. 
done

RUN=ENDOPHYTE
for s in BAC FUN; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mv $PROJECT_FOLDER/data/$s/filtered/*E* $PROJECT_FOLDER/data/$RUN/$s/filtered/. 
  mv $PROJECT_FOLDER/data/$s/unfiltered/*E* $PROJECT_FOLDER/data/$RUN/$s/unfiltered/. 
done

