# folder containing project files
PROJECT_FOLDER=~/projects/Endophytes

# folder for specific sequencer run
RUN=2016

# make the project folder
mkdir -p $PROJECT_FOLDER

# link to metabarcoding pipeline (specific for authors profile)
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

# folder to hold fatsq files
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq

# variable to hold folder names (BAC and FUN)
RIB="BAC FUN"

# loop through the RIB variable, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for s in $RIB; do
mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# This data is not multiplexed so fastq files can be decompressed to BAC and FUN folders 

for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*.gz; do 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c unzip $FILE 2
done

cd $PROJECT_FOLDER/data/$RUN/fastq

# move files into correct place
mv Apple-Lady* ../BAC/fastq/.
mv Apple-* ../FUN/fastq/.

# oops move gz files back again
mv ../BAC/fastq/*.gz .
mv ../FUN/fastq/*.gz .

# pre-process BAC files (min length 300, min overlap 5, quality 0.5)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/BAC \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 5 0.5

# Pre-process FUN files (min length 200, MAX R2 length 250, quality 1)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 200 250 1

# identify none ITS regions (R1 only)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
 $PROJECT_FOLDER/data/$RUN/FUN/fasta \
 R1 \
 $PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
 $PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
 ssu 58ss 20

# remove none ITS regions from sequence
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$PROJECT_FOLDER/data/$RUN/FUN/fasta/*R1" \
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58"

# Move merged fasta to filtered folder
for f in $PROJECT_FOLDER/data/$RUN/$SSU/unfiltered/*r1*; do 
 F=$(echo $f|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $f|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") ; 
 mv ../fasta/${L}_R1/$F ../filtered/$L; done


### 2017 DATA ###

# folder for specific sequencer run
RUN=2017

# loop through the RIB variable, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for s in $RIB; do
mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# data was already pre-processed as part of Kenyan beans project
# manual copy of BAC/FUN data to filtered/unfiltered folders


### COMBINED DATA ###

RUN=COMBINED

# loop through the RIB variable, i.e. s = BAC on first loop, S= FUN on second loop, and create the folders
for s in $RIB; do
mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# Symlink files to COMBINED folders
for s in $RIB; do
ln -s $PROJECT_FOLDER/data/2016/$s/filtered/* $PROJECT_FOLDER/data/$RUN/$s/filtered/.
ln -s $PROJECT_FOLDER/data/2017/$s/filtered/* $PROJECT_FOLDER/data/$RUN/$s/filtered/.
ln -s $PROJECT_FOLDER/data/2016/$s/unfiltered/* $PROJECT_FOLDER/data/$RUN/$s/unfiltered/.
ln -s $PROJECT_FOLDER/data/2017/$s/unfiltered/* $PROJECT_FOLDER/data/$RUN/$s/unfiltered/.
done

