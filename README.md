# apple_replant
Metagenomic study of apple replant disease
####Installing Qiime to a local directory
Downloaded Python 2.7.9 tar ball and unzipped
From root of Python 2.7.9 directory ran :

	./configure --prefix=/home/deakig/usr/local --exec-prefix=/home/deakig/usr/local --enable-unicode=ucs4
	make
	make install

Downloaded pip tarball amd unzipped to pip directory then ran

	~/usr/local/bin/python ~/pip/getpip.py


Set qiime path with (not permanent)

	export PYTHONUSERBASE=/home/deakig/usr/local/
	
To install Qiime and dependencies

	~/usr/local/bin/python -m pip install --user --upgrade --force-reinstall numpy
	~/usr/local/bin/python -m pip install --user --upgrade --force-reinstall qiime
	
(the upgrade and force-reinstall flags are probably not necessary)

To test qiime, ensure ~/usr/local/bin (the qiime script directory) is in path

	export PATH=$PATH:/home/deakig/usr/local/bin

then

	 ~/usr/local/bin/python ~/usr/local/bin/print_qiime_config.py -t

should retun something like 
$>Ran 9 test in 0.05s
$>OK

###Parallel qiime
for single machine throw in -a -O (no. processes) to the workflow script

using HPC...

	create qimme_config in home root
	cd ~
	touch .qiime_config

added:
jobs_to_start 8
temp_dir $HOME/tmp
cluster_jobs_fp start_parallel_jobs_sc.py	

hacked start_parallel_jobs_sc.py for use in our environment

#16s workflow

###Join PE reads
Change to trimmed directory then run below script (this will also do ITS samples)

	for f in ./*trimmed*; 
		do counter=$((counter+1)); 
		if (( $counter % 2 == 0 )); 
			then R2=$f;
			echo join_paired_ends.py -f $R1 -r $R2 -o $counter;
			join_paired_ends.py -f $R1 -r $R2 -o $counter; 
		fi; 
		R1=$f; 
	done


###rename files 

	counter=84
	for d in * 
		do counter=$((counter+1));
		cd S$counter
		for f in *
		do
			mv -i "${f}" "S${f/fastqjoin/$counter}"
		done
		cd ..
	done

###convert joined fastq to fasta

	counter=84
	for d in *
		do counter=$((counter+1));
		cd S$counter
		for f in ./*join*
		do
			../../../scripts/fq2fa.pl $f $f.fa S$d
			mv $f.fa ../../fasta/.
		done
		cd ..
	done

####concatenate files
created 16S and ITS under fasta folder and moved file to appropriate place

	cat S* >16S.joined.fa	

####OTU Picking and descriptive statistics

	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/16S/16S.joined.fa /home/deakig/projects/metagenomics/analysis/16S_otus /home/deakig/projects/metagenomics/scripts/parameters.txt /home/deakig/usr/local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta TRUE

	 X=`biom summarize-table -i analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"`

	./core_diversity.sh /home/deakig/projects/metagenomics/analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom /home/deakig/projects/metagenomics/analysis/16s_cdout/ /home/deakig/projects/metagenomics/data/map.tsv /home/deakig/projects/metagenomics/analysis/16S_otus/rep_set.tre $X

#### Statistical analysis
	
	Rscript analysis.R "analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom" 6 median res.sig.csv


#OLD Stuff 
###trim trimmomatic
```shell
./trimmomatic.sh Replant-1A_S14_L001_R1_001.fastq Replant-1A_S14_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-1A_S30_L001_R1_001.fastq Replant-1A_S30_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-5A_S15_L001_R1_001.fastq Replant-5A_S15_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-1A_S38_L001_R1_001.fastq Replant-1A_S38_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-5A_S39_L001_R1_001.fastq Replant-5A_S39_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-1A_S22_L001_R1_001.fastq Replant-1A_S22_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-5A_S23_L001_R1_001.fastq Replant-5A_S23_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-1A_S6_L001_R1_001.fastq Replant-1A_S6_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-5A_S7_L001_R1_001.fastq Replant-5A_S7_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
./trimmomatic.sh Replant-5A_S31_L001_R1_001.fastq Replant-5A_S31_L001_R2_001.fastq /home/deakig/projects/metagenomics/data /home/deakig/projects/metagenomics/scripts
```
##join paired end reads (pretty ugly as it bungs them into a unique folder but files have same name...
```shell
join_paired_ends.py -f Replant-1A_S14_L001_R1_001.fastq.trimmed.fq -r Replant-1A_S14_L001_R2_001.fastq.trimmed.fq -o Replant-1A_S14
join_paired_ends.py -f Replant-1A_S30_L001_R1_001.fastq.trimmed.fq -r Replant-1A_S30_L001_R2_001.fastq.trimmed.fq -o Replant-1A_S30
join_paired_ends.py -f Replant-5A_S15_L001_R1_001.fastq.trimmed.fq -r Replant-5A_S15_L001_R2_001.fastq.trimmed.fq -o Replant-5A_S15
join_paired_ends.py -f Replant-1A_S38_L001_R1_001.fastq.trimmed.fq -r Replant-1A_S38_L001_R2_001.fastq.trimmed.fq -o Replant-1A_S38
join_paired_ends.py -f Replant-5A_S39_L001_R1_001.fastq.trimmed.fq -r Replant-5A_S39_L001_R2_001.fastq.trimmed.fq -o Replant-5A_S39
join_paired_ends.py -f Replant-1A_S22_L001_R1_001.fastq.trimmed.fq -r Replant-1A_S22_L001_R2_001.fastq.trimmed.fq -o Replant-1A_S22
join_paired_ends.py -f Replant-5A_S23_L001_R1_001.fastq.trimmed.fq -r Replant-5A_S23_L001_R2_001.fastq.trimmed.fq -o Replant-5A_S23
join_paired_ends.py -f Replant-1A_S6_L001_R1_001.fastq.trimmed.fq -r Replant-1A_S6_L001_R2_001.fastq.trimmed.fq -o Replant-1A_S6
join_paired_ends.py -f Replant-5A_S7_L001_R1_001.fastq.trimmed.fq -r Replant-5A_S7_L001_R2_001.fastq.trimmed.fq -o Replant-5A_S7
join_paired_ends.py -f Replant-5A_S31_L001_R1_001.fastq.trimmed.fq -r Replant-5A_S31_L001_R2_001.fastq.trimmed.fq -o Replant-5A_S31
```
####cat files (probably best)
	cat  Replant-1A_S14/* > Replant-1A_S14.all.fq
	cat  Replant-1A_S30/* > Replant-1A_S30.all.fq
	cat  Replant-5A_S15/* > Replant-5A_S15.all.fq
	cat  Replant-1A_S38/* > Replant-1A_S38.all.fq
	cat  Replant-5A_S39/* > Replant-5A_S39.all.fq
	cat  Replant-1A_S22/* > Replant-1A_S22.all.fq
	cat  Replant-5A_S23/* > Replant-5A_S23.all.fq
	cat  Replant-1A_S6/* > Replant-1A_S6.all.fq
	cat  Replant-5A_S7/* > Replant-5A_S7.all.fq
	cat  Replant-5A_S31/* > Replant-5A_S31.all.fq



####convert to fasta
	../../scripts/fq2fa.pl Replant-1A_S14.all.fq Replant-1A_S14.all.fa Replant-1A_S14
	../../scripts/fq2fa.pl Replant-1A_S30.all.fq Replant-1A_S30.all.fa Replant-1A_S30
	../../scripts/fq2fa.pl Replant-5A_S15.all.fq Replant-5A_S15.all.fa Replant-5A_S15
	../../scripts/fq2fa.pl Replant-1A_S38.all.fq Replant-1A_S38.all.fa Replant-1A_S38
	../../scripts/fq2fa.pl Replant-5A_S39.all.fq Replant-5A_S39.all.fa Replant-5A_S39
	../../scripts/fq2fa.pl Replant-1A_S22.all.fq Replant-1A_S22.all.fa Replant-1A_S22
	../../scripts/fq2fa.pl Replant-5A_S23.all.fq Replant-5A_S23.all.fa Replant-5A_S23
	../../scripts/fq2fa.pl Replant-1A_S6.all.fq Replant-1A_S6.all.fa Replant-1A_S6
	../../scripts/fq2fa.pl Replant-5A_S7.all.fq Replant-5A_S7.all.fa Replant-5A_S7
	../../scripts/fq2fa.pl Replant-5A_S31.all.fq Replant-5A_S31.all.fa Replant-5A_S31

---- update
	cat * > fasta/all.fa

###create  qiime mapping file
(hash must be present before SampleID) 
#SampleID	BarcodeSequence	LinkerPrimerSequence	SampleType	Description
Replant.1A.S14.L001.R1			1	type 1
Replant.1A.S30.L001.R2			1	type 1
Replant.5A.S15.L001.R1			5	type 5
Replant.5A.S31.L001.R2			5	type 5
Replant.1A.S14.L001.R2			1	type 1
Replant.1A.S38.L001.R1			1	type 1
Replant.5A.S15.L001.R2			5	type 5
Replant.5A.S39.L001.R1			5	type 5
Replant.1A.S22.L001.R1			1	type 1
Replant.1A.S38.L001.R2			1	type 1
Replant.5A.S23.L001.R1			5	type 5
Replant.5A.S39.L001.R2			5	type 5
Replant.1A.S22.L001.R2			1	type 1
Replant.1A.S6.L001.R1			1	type 1
Replant.5A.S23.L001.R2			5	type 5
Replant.5A.S7.L001.R1			5	type 5
Replant.1A.S30.L001.R1			1	type 1
Replant.1A.S6.L001.R2			1	type 1
Replant.5A.S31.L001.R1			5	type 5
Replant.5A.S7.L001.R2			5	type 5

###create qiime parameters file (this will need to be modified depending on the data type and analyses required)


###quiime commands

#OTU picker
	pick_open_reference_otus.py -f -o otus -i data/fasta/all.fa -p scripts/parameters.txt 

#summerise data (note min sequencing depth - or drop samples)
	biom summarize-table -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom

#diversity analysis
	core_diversity_analyses.py -o cdout/ -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom -m data/map.tsv -t otus/rep_set.tre -e 71941 --suppress_beta_diversity
#-e is sequencing depth from summarise table, --suppress_beta_diversity suppresses the emperor 3d pca plots (useful if less than 4 samples).



	./core_diversity.sh /home/deakig/projects/metagenomics/otus2/otu_table_mc2_w_tax_no_pynast_failures.biom /home/deakig/projects/metagenomics/cdout_cluster/ /home/deakig/projects/metagenomics/data/map2.tsv /home/deakig/projects/metagenomics/otus2/rep_set.tre 3272
#to do - write cluster version of OTU picker
done...
	./pick_OTU.sh DATA OUT_FOLDER PARAM_FILE REF_FILE
default ref file is something like:
	/home/deakig/usr/local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta


###non-default reference database
#add to parameters file:
assign_taxonomy:id_to_taxonomy_fp /path_to_taxonomy_txt
assign_taxonomy:reference_seqs_fp /path_to_fasta

	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/all.fa /home/deakig/projects/metagenomics/debug/otus_16_18_taxa /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/Silva119/97/Silva_119_rep_set97_aligned.fna
biom summarize-table -i debug/otus_16_18_taxa/otu_table_mc2_w_tax_no_pynast_failures.biom

	./core_diversity.sh /home/deakig/projects/metagenomics/debug/otus_16_18_taxa/otu_table_mc2_w_tax_no_pynast_failures.biom /home/deakig/projects/metagenomics/debug/cdout_16_18_taxa /home/deakig/projects/metagenomics/data/map2.tsv /home/deakig/projects/metagenomics/debug/otus_16_18_taxa/rep_set.tre 3357


/home/deakig/projects/metagenomics/taxonomies/Silva119/97/Silva_119_rep_set97_aligned.fna


########NEW DATA#########
#changed trimmomatic script to use all files in directory
	./trimmomatic.sh /home/deakig/projects/metagenomics/data/replant2 /home/deakig/projects/metagenomics/scripts

	for f in ./*trimmed*; 
	do counter=$((counter+1)); 
		if (( $counter % 2 == 0 )); 
			then R2=$f;
			echo join_paired_ends.py -f $R1 -r $R2 -o $counter;
			join_paired_ends.py -f $R1 -r $R2 -o $counter; 
		fi; 
	R1=$f; 
	done

	cat S85/* > S85.all.fq
	cat S86/* > S86.all.fq
	cat S87/* > S87.all.fq
	cat S88/* > S88.all.fq
	cat S89/* > S89.all.fq
	cat S90/* > S90.all.fq
	cat S91/* > S91.all.fq
	cat S92/* > S92.all.fq
	cat S93/* > S93.all.fq
	cat S94/* > S94.all.fq
	cat S95/* > S95.all.fq
	cat S96/* > S96.all.fq
	
	../../scripts/fq2fa.pl S85.all.fq S85.all.fa S85
	../../scripts/fq2fa.pl S86.all.fq S86.all.fa S86
	../../scripts/fq2fa.pl S87.all.fq S87.all.fa S87
	../../scripts/fq2fa.pl S88.all.fq S88.all.fa S88
	../../scripts/fq2fa.pl S89.all.fq S89.all.fa S89
	../../scripts/fq2fa.pl S90.all.fq S90.all.fa S90
	../../scripts/fq2fa.pl S91.all.fq S91.all.fa S91
	../../scripts/fq2fa.pl S92.all.fq S92.all.fa S92
	../../scripts/fq2fa.pl S93.all.fq S93.all.fa S93
	../../scripts/fq2fa.pl S94.all.fq S94.all.fa S94
	../../scripts/fq2fa.pl S95.all.fq S95.all.fa S95
	../../scripts/fq2fa.pl S96.all.fq S96.all.fa S96
	
	cat S* > all_r2.fa
###16S
./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/all_r2.fa /home/deakig/projects/metagenomics/analysis/otus /home/deakig/projects/metagenomics/scripts/parameters.txt /home/deakig/usr/local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta

biom summarize-table -i analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom

./core_diversity.sh /home/deakig/projects/metagenomics/analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom /home/deakig/projects/metagenomics/analysis/cdout/ /home/deakig/projects/metagenomics/data/map.tsv /home/deakig/projects/metagenomics/analysis/otus/rep_set.tre 197641


./core_diversity.sh /home/deakig/projects/metagenomics/analysis/otus/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/analysis/test_16s/ /home/deakig/projects/metagenomics/data/map.tsv  . 137672

###Fungal
./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/all_r2.fa /home/deakig/projects/metagenomics/analysis/otus /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta


####Trimming testing with bacterial and fungal primers
S85 and S91
default timming: 
S85	joined 		145740
   	unjoined	89642
S91 	joined		15267
    	unjoined	156266	

default + headcrop:9
S85	joined		145640
	unjoined	89574
S91	joined		12278
	unjoined	158639

default + headcrop:9 + maxlength 140
rubbish

Using default
######

######Testing Fungal for F, R, PE, PE + F_un, PE + R_un, PE + F_un + R_un
#convert F to fa
../../scripts/fq2fa.pl 91_S91_L001_R1_001.fastq.gz.trimmed.fq S91.f.fa S91
../../scripts/fq2fa.pl 92_S92_L001_R1_001.fastq.gz.trimmed.fq S92.f.fa S92
../../scripts/fq2fa.pl 93_S93_L001_R1_001.fastq.gz.trimmed.fq S93.f.fa S93
../../scripts/fq2fa.pl 94_S94_L001_R1_001.fastq.gz.trimmed.fq S94.f.fa S94
../../scripts/fq2fa.pl 95_S95_L001_R1_001.fastq.gz.trimmed.fq S95.f.fa S95
../../scripts/fq2fa.pl 96_S96_L001_R1_001.fastq.gz.trimmed.fq S96.f.fa S96

#convert R to fa
../../scripts/fq2fa.pl 91_S91_L001_R2_001.fastq.gz.trimmed.fq S91.r.fa S91
../../scripts/fq2fa.pl 92_S92_L001_R2_001.fastq.gz.trimmed.fq S92.r.fa S92
../../scripts/fq2fa.pl 93_S93_L001_R2_001.fastq.gz.trimmed.fq S93.r.fa S93
../../scripts/fq2fa.pl 94_S94_L001_R2_001.fastq.gz.trimmed.fq S94.r.fa S94
../../scripts/fq2fa.pl 95_S95_L001_R2_001.fastq.gz.trimmed.fq S95.r.fa S95
../../scripts/fq2fa.pl 96_S96_L001_R2_001.fastq.gz.trimmed.fq S96.r.fa S96


##F
#cat files
cat *.f.fa > fungal.f.fa
#R
#cat files
cat *.r.fa > fungal.r.fa
#PE
#cat files
cat *.PE.fa > fungal.PE.fa
#PE + F
#cat files
cat *.pef.fa > fungal.pef.fa
#PE + R
#cat files
cat *.per.fa > fungal.per.fa
#PE + F + R
#cat files
cat *.all.fa > fungal.all.fa

##OTU picking
../scripts/pick_OTU.sh  /home/deakig/projects/metagenomics/testing/fungal.f.fa /home/deakig/projects/metagenomics/testing/f_f /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta
../scripts/pick_OTU.sh  /home/deakig/projects/metagenomics/testing/fungal.r.fa /home/deakig/projects/metagenomics/testing/f_r /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta
../scripts/pick_OTU.sh  /home/deakig/projects/metagenomics/testing/fungal.PE.fa /home/deakig/projects/metagenomics/testing/f_PE /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta
../scripts/pick_OTU.sh  /home/deakig/projects/metagenomics/testing/fungal.pef.fa /home/deakig/projects/metagenomics/testing/f_pef /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta
../scripts/pick_OTU.sh  /home/deakig/projects/metagenomics/testing/fungal.per.fa /home/deakig/projects/metagenomics/testing/f_per /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta
../scripts/pick_OTU.sh  /home/deakig/projects/metagenomics/testing/fungal.all.fa /home/deakig/projects/metagenomics/testing/f_all /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta

#basic stats
biom summarize-table -i f_x/otu_table_mc2_w_tax.biom

../scripts/core_diversity.sh /home/deakig/projects/metagenomics/testing/f_f/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/testing/cdout_f_f/ /home/deakig/projects/metagenomics/data/map.tsv  . 75961
../scripts/core_diversity.sh /home/deakig/projects/metagenomics/testing/f_r/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/testing/cdout_f_r/ /home/deakig/projects/metagenomics/data/map.tsv  . 64712
../scripts/core_diversity.sh /home/deakig/projects/metagenomics/testing/f_PE/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/testing/cdout_f_PE/ /home/deakig/projects/metagenomics/data/map.tsv  . 2876
../scripts/core_diversity.sh /home/deakig/projects/metagenomics/testing/f_pef/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/testing/cdout_f_pef/ /home/deakig/projects/metagenomics/data/map.tsv  . 75860
../scripts/core_diversity.sh /home/deakig/projects/metagenomics/testing/f_per/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/testing/cdout_f_per/ /home/deakig/projects/metagenomics/data/map.tsv  . 64912
../scripts/core_diversity.sh /home/deakig/projects/metagenomics/testing/f_all/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/testing/cdout_f_all/ /home/deakig/projects/metagenomics/data/map.tsv  . 137925


###Paired-end updates
May be better to junk reads less than about 150 bp


###Removal of ssu and 5.8 from ITS
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S91/fastqjoin.un1.fastq S91.un.1.fa S91
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S91/fastqjoin.un2.fastq S91.un.2.fa S91
~projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S91/fastqjoin. S91.un.2.fa S91
fastqjoin.join.fastq  fastqjoin.un1.fastq   fastqjoin.un2.fastq   
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S91/fastqjoin.join.fastq S91.pe.fa S91PE
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S92/fastqjoin.un1.fastq S92.un.1.fa S92
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S93/fastqjoin.un1.fastq S93.un.1.fa S93
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S94/fastqjoin.un1.fastq S94.un.1.fa S94
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S95/fastqjoin.un1.fastq S95.un.1.fa S95
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S96/fastqjoin.un1.fastq S96.un.1.fa S96
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S92/fastqjoin.un2.fastq S92.un.2.fa S92
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S93/fastqjoin.un2.fastq S93.un.2.fa S93
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S94/fastqjoin.un2.fastq S94.un.2.fa S94
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S95/fastqjoin.un2.fastq S95.un.2.fa S95
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S96/fastqjoin.un2.fastq S96.un.2.fa S96
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S91/fastqjoin.join.fastq S91.pe.fa S91
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S92/fastqjoin.join.fastq S92.pe.fa S92
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S93/fastqjoin.join.fastq S93.pe.fa S93
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S94/fastqjoin.join.fastq S94.pe.fa S94
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S95/fastqjoin.join.fastq S95.pe.fa S95
~/projects/metagenomics/data/joined$ ../../scripts/fq2fa.pl S96/fastqjoin.join.fastq S96.pe.fa S96
