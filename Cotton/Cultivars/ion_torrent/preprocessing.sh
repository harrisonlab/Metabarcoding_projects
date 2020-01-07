PROJECT_FOLDER=~/projects/Cotton/
mkdir -p $PROJECT_FOLDER
ln -s $MBPL $PROJECT_FOLDER/metabarcoding_pipeline 

RUN=cultivar_ion

for s in BAC FUN; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done


# QC
for FILE in $PROJECT_FOLDER/data/$RUN/BAC/fastq/*; do
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/BAC/quality
done

for FILE in $PROJECT_FOLDER/data/$RUN/FUN/fastq/*; do
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/FUN/quality
done

# Length trimming
#cat $PROJECT_FOLDER/$RUN/BAC/fastq/* > bac.cat.fq
#cat $PROJECT_FOLDER/$RUN/FUN/fastq/* > fun.cat.fq

## ion torrent S5 looks like it uses extra phred 33 charcters (L and M - maybe more?) below to check
# awk 'NR % 4 ==0' LM28.D10.fastq|tr -d '\n'|grep -o . |sort -u|paste -s -d '\0'
# cat xaa|tr LM K > xaa1


#bacteria
SSU=BAC
QUAL=0.005
MAX_LEN=400

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ion \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
$QUAL $MAX_LEN

#Fungi
SSU=FUN
QUAL=0.01
MAX_LEN=150

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ion \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
$QUAL $MAX_LEN

