# Create folders
PROJECT_FOLDER=~/projects/ARD
RUN=160711
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous

RIB="BAC FUN"
for s in $RIB; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done

# demultiplex 
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

# Bacteria workflow
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
 "$PROJECT_FOLDER/data/$RUN/BAC/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/BAC \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 5 0.5 17 21 

# Fungi workflow
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 200 1 23 21

for F in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*_R1.fa; do 
  FO=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
  L=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}') ;
  awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$FO;
done

# Create Koekoek folders and links 
RIB="BAC FUN"
for s in $RIB; do
  mkdir -p $PROJECT_FOLDER/analysis/Koekoek/$s/filtered
  ln -s $PROJECT_FOLDER/data/$RUN/$s/filtered/K-* $PROJECT_FOLDER/analysis/Koekoek/$s/filtered/.
  mkdir -p $PROJECT_FOLDER/analysis/Koekoek/$s/unfiltered
  ln -s $PROJECT_FOLDER/data/$RUN/$s/unfiltered/K-* $PROJECT_FOLDER/analysis/Koekoek/$s/unfiltered/.
done
