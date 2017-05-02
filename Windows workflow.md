## Introduction
Because I've been asked...

## Demultiplexing
```
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
"$ARDERI/data/$RUN/fastq/*_R1_*" 0 \
$P1F $P1R $P2F $P2R
```
## Preprocessing


## ITS specific (if wanted)

## UPARSE pipeline

### Cluster
```
#denoise
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $ARDERI $RUN $SSU $FPL $RPL

# or clustering with cluser_otu
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UCLUS $ARDERI $RUN $SSU $FPL $RP
```

### Assign taxonomy
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```

### Create OTU table
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```
