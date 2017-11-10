# setup links and variables
PROJECT_FOLDER=~/projects/ARD/analysis
ln -s ~/pipelines/metabarcoding/ $PROJECT_FOLDER/metabarcoding_pipeline
RUN="../Rootstocks/Rhizosphere"
mkdir ~/projects/ARD/analysis/data # I've hardcoded some stuff in the pipeline to use a folder called data (need to edit this...) 

# Fungi
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/FUN/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/FUN \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 200 1 23 21

for F in $PROJECT_FOLDER/data/$RUN/FUN/fasta/*_R1.fa; do 
 FO=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}') ;
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/FUN/filtered/$FO;
done

# Bacteria
SSU=BAC;FPL=17;RPL=21;MINL=300;MINOVER=5;QUAL=0.5
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
$PROJECT_FOLDER/data/$RUN/$SSU \
$PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL $FPL $RPL 
