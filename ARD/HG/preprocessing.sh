# sequence data location(s)
RUN=170224
RUN=170323

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

# bacteria
SSU=BAC;FPL=17;RPL=21
MINL=300;MINOVER=5;QUAL=0.5

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
$PROJECT_FOLDER/data/$RUN/$SSU \
$PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL

# fungi
SSU=FUN;FPL=23;RPL=21
MINL=200;MAXR2=250;QUAL=1

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 $MINL $MAXR2 $QUAL; 

## forward only pipeline to remove none ITS regions
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
$PROJECT_FOLDER/data/$RUN/$SSU/fasta \
R1 \
$PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
$PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
ssu 58ss 20

### 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fasta/*R1" \
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58"
 
###
for D in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*1; do 
 F=$(echo $D|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $D|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") ;
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $D/$F > $F.tmp && mv $F.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$F;
# mv $PROJECT_FOLDER/data/$RUN/$SSU/fasta/${L}_R1/$F $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$L; 
done
