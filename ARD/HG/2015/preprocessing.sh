# folder containing project files
PROJECT_FOLDER=~/projects/ARD

# sequencer run folders 
RUNS="151019 160408 160418 160506 160111 160122 160201" 

# make the project folder
mkdir -p $PROJECT_FOLDER

# link to metabarcoding pipeline (specific for authors profile)
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

# folder to hold fatsq files
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq

# variable to hold folder names (BAC and FUN)
RIB="16S ITS"

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
   300 5 0.5
done

# Pre-process ITS files (min length 200, MAX length 300, quality 1)
for RUN in $RUNS; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
   "$PROJECT_FOLDER/data/$RUN/ITS/fastq/*R1*.fastq" \
   $PROJECT_FOLDER/data/$RUN/ITS \
   $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
   200 300 1
done

## identify none ITS (FUN) regions
for RUN in $RUNS; do
  # forward reads
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
   $PROJECT_FOLDER/data/$RUN/ITS/fasta \
   R1 \
   $PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
   $PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
   ssu 58ss 20
  # reverse reads 
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
   $PROJECT_FOLDER/data/$RUN/ITS/fasta \
   R2 \
   $PROJECT_FOLDER/metabarcoding_pipeline/hmm/lsu_start.hmm \
   $PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_end.hmm \
   lsu    
done

## remove none ITS regions from sequence
for RUN in $RUNS; do
  # forward reads
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
   "$PROJECT_FOLDER/data/$RUN/ITS/fasta/*R1" \
   $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
   "*.\\.ssu" \
   "*.\\.58"
  # reverse reads
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
   "$PROJECT_FOLDER/data/$RUN/ITS/fasta/*R2" \
   $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_58Se_LSU_v2.R \
   "*.\\.58" \
   "*.\\.lsu" \
   FALSE   
done

# merge forward and reverse reads
for RUN in $RUNS; do
  # move files to filtered folder
  find $PROJECT_FOLDER/data/$RUN/ITS/fasta -type f -name *.r*.fa|xargs -I myfile mv myfile $PROJECT_FOLDER/data/$RUN/ITS/filtered/.
  cd $PROJECT_FOLDER/data/$RUN/$SSU/filtered
  for f in $PROJECT_FOLDER/data/$RUN/$SSU/filtered/*r1.fa
  do
      R1=$f
      R2=$(echo $R1|sed 's/\.r1\.fa/\.r2\.fa/')
      S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
      $PROJECT_FOLDER/metabarcoding_pipeline/scripts/catfiles_v2.pl $R1 $R2 $S;
  done
  mkdir R1
  mkdir R2
  mv *r1* R1/.
  mv *r2* R2/.
done
