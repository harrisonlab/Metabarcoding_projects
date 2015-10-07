# apple_replant
Metagenomic study of apple replant disease
## Installing Qiime to a local directory
Downloaded Python 2.7.9 tar ball and unzipped.

From root of Python 2.7.9 directory ran :

	./configure --prefix=/home/deakig/usr/local --exec-prefix=/home/deakig/usr/local --enable-unicode=ucs4
	make
	make install

Downloaded pip tarball amd unzipped to pip directory then ran:

	~/usr/local/bin/python ~/pip/getpip.py


Set Qiime path with (not permanent)

	export PYTHONUSERBASE=/home/deakig/usr/local/
	
To install Qiime and dependencies

	~/usr/local/bin/python -m pip install --user --upgrade --force-reinstall numpy
	~/usr/local/bin/python -m pip install --user --upgrade --force-reinstall qiime
	
(the upgrade and force-reinstall flags may not be necessary)

To test qiime, ensure ~/usr/local/bin (the qiime script directory) is in path

	export PATH=$PATH:/home/deakig/usr/local/bin

then

	 ~/usr/local/bin/python ~/usr/local/bin/print_qiime_config.py -t

should retun something like

	$> Ran 9 test in 0.05s
	$> OK

### Parallel qiime
for single machine throw in -a -O (no. processes) to the workflow script

using HPC...

	create qimme_config in home root
	cd ~
	touch .qiime_config

added:
jobs_to_start 8
temp_dir $HOME/tmp
cluster_jobs_fp start_parallel_jobs_sc.py	

hacked start_parallel_jobs_sc.py for use in our environment.
___
## Common workflow
### QC
Qualtiy checking was performed with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

From same folder containing fastq files ran:

	fastqc *

### Trimming
Paired end trimming was preformed with Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic).

The following settings were used:
- minimum length 200
- sliding window
 + 8 bases
 + quality 15
- illumina adapter clipping
 + 2 mismatches 
 + palindrome quality 30
 + clip threshold quality 10

Shell script trim.sh used to submit trimming jobs to cluster. 

The first argument specifies a directory of paired end fastq files (all file in folder will be processed - paired files must sort adjacently)

Second argument specifies location of illumina adapter file
```shell	
./trim.sh /home/deakig/projects/metagenomics/data/fastq /home/deakig/projects/metagenomics/scripts 
```
## 16s workflow

### Join PE reads
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


### Rename files 
The counter used in the next couple of commands was set to match the names of the samples, i.e. S85, S86 and etc.

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

### Convert joined fastq to fasta

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

### Remove chimeras
downloaded usearch 8.0 and RDP gold reference database from http://drive5.com/usearch/manual/cmd_uchime_ref.html

	counter=84
	for f in /home/deakig/projects/metagenomics/data/fasta/16S/*
	do counter=$((counter+1));
		./chimeras.sh $f /home/deakig/projects/metagenomics/taxonomies/RDP_gold.fasta S${counter}.cfree.fa
		/home/deakig/projects/metagenomics/data/fasta/de_chimeraed/
	done

#### Concatenate files
	cat *cfree* >16S.joined.fa	

### OTU Picking and descriptive statistics

	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/16S/16S.joined.fa /home/deakig/projects/metagenomics/analysis/16S_otus /home/deakig/projects/metagenomics/scripts/parameters.txt /home/deakig/usr/local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta TRUE
	 X=`biom summarize-table -i analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"`
	./core_diversity.sh /home/deakig/projects/metagenomics/analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom /home/deakig/projects/metagenomics/analysis/16s_cdout/ /home/deakig/projects/metagenomics/data/map.tsv /home/deakig/projects/metagenomics/analysis/16S_otus/rep_set.tre $X

### Statistical analysis
Requires a colData file describing condition (e.g. infected or uninfected) for each sample
analysis.R biom_table "no. samples" median/geomean outfile
	
	Rscript analysis.R "analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom" 6 median res.sig.csv
	
## ITS workflow

### Convert to unpaired fasta files

	X=91
	counter=0
	for f in ./*trimmed*;
	do counter=$((counter+1));
	  if [ "$counter" -gt 12 ]
	  then
	    if (( $counter % 2 == 0 ))
	    then
	      ../../scripts/fq2fa.pl $f ../fasta/ITS/${f}.fa S$X ;
	      X=$((X+1))
	    else
	      ../../scripts/fq2fa.pl $f ../fasta/ITS/${f}.fa S$X ;
	    fi
	  fi
	done

### Rename files 

	X=91
	counter=0
	mkdir -p new
	for f in *trimmed*;
	do counter=$((counter+1));
	  if (( $counter % 2 == 0 ))
	  then
	    mv -i "${f}" "S${X}_R2.fa"
	    X=$((X+1))
	  else
	    mv -i "${f}" "S${X}_R1.fa"
	  fi
	done

### SSU/58S/LSU removal
Using HHMMER v 3.x from ....

Download HMM files from ITSx (need website)

Hacked the HMM files to include a MAXL satement (required) and split out SSU,58S and LSU into seperate files (using fungal only)

	perl cut_hmm v.3.1 chopped_hmm fungi
	cd chopped_hmm
	cat *SSU*> ssu_end.hmm
	cat *58S_start* > 58s_start.hmm
	cat *58S_end* > 58s_end.hmm
	cat *LSU* > lsu_start.hmm
	hmmpress ssu_end.hmm
	hmmpress 58s_end.hmm
	hmmpress 58s_start.hmm
	hmmpress lsu_start.hmm
	
##### Split file into chunks for SSU/58S/LSu removal
	X=91 
	counter=0
	for f in *.fa;
	do counter=$((counter+1));
	  if (( $counter % 2 == 0 ))
	  then
	    mkdir S${X}_R2
	    split -l 2000 $f -a 3 -d S${X}_R2/$f.
	    X=$((X+1))
	  else
	    mkdir S${X}_R1
	    split -l 2000 $f -a 3 -d S${X}_R1/$f.
	  fi
	done

##### Remove SSU/LSU

	counter=0
	for d in */;
	do counter=$((counter+1));
	  cd $d
	  if (( $counter % 2 == 0 ))
	  then
	    for f in *;
	    do 
	      nscan.sh /home/deakig/projects/metagenomics/data/fasta/ITS/${d}$f /home/deakig/projects/metagenomics/data/fasta/ITS/${d}${f}.lsu 20 /home/deakig/projects/metagenomics/hmm/lsu_start.hmm
	      nscan.sh /home/deakig/projects/metagenomics/data/fasta/ITS/${d}$f /home/deakig/projects/metagenomics/data/fasta/ITS/${d}${f}.58se 20 /home/deakig/projects/metagenomics/hmm/58s_end.hmm
	    done
	  else
	    for f in *;
	    do
	      nscan.sh /home/deakig/projects/metagenomics/data/fasta/ITS/${d}$f /home/deakig/projects/metagenomics/data/fasta/ITS/${d}${f}.ssu 20 /home/deakig/projects/metagenomics/hmm/ssu_end.hmm
	      nscan.sh /home/deakig/projects/metagenomics/data/fasta/ITS/${d}$f /home/deakig/projects/metagenomics/data/fasta/ITS/${d}${f}.58ss 20 /home/deakig/projects/metagenomics/hmm/58s_start.hmm
	    done
	  fi
	  cd ..
	done

##### Merge output
	./ITS.sh /home/deakig/projects/metagenomics/rm_SSU_58Ss.R /home/deakig/projects/metagenomics/data/fasta/ITS/S91_R1/ "*.\\.ssu" "*.\\.58" /home/deakig/projects/metagenomics/data/fasta/ITS/S91_R1.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_SSU_58Ss.R /home/deakig/projects/metagenomics/data/fasta/ITS/S92_R1/ "*.\\.ssu" "*.\\.58" /home/deakig/projects/metagenomics/data/fasta/ITS/S92_R1.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_SSU_58Ss.R /home/deakig/projects/metagenomics/data/fasta/ITS/S93_R1/ "*.\\.ssu" "*.\\.58" /home/deakig/projects/metagenomics/data/fasta/ITS/S93_R1.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_SSU_58Ss.R /home/deakig/projects/metagenomics/data/fasta/ITS/S94_R1/ "*.\\.ssu" "*.\\.58" /home/deakig/projects/metagenomics/data/fasta/ITS/S94_R1.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_SSU_58Ss.R /home/deakig/projects/metagenomics/data/fasta/ITS/S95_R1/ "*.\\.ssu" "*.\\.58" /home/deakig/projects/metagenomics/data/fasta/ITS/S95_R1.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_SSU_58Ss.R /home/deakig/projects/metagenomics/data/fasta/ITS/S96_R1/ "*.\\.ssu" "*.\\.58" /home/deakig/projects/metagenomics/data/fasta/ITS/S96_R1.fa
	
	./ITS.sh /home/deakig/projects/metagenomics/rm_58Se_LSU.R /home/deakig/projects/metagenomics/data/fasta/ITS/S91_R2/ "*.\\.58" "*.\\.lsu" /home/deakig/projects/metagenomics/data/fasta/ITS/S91_R2.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_58Se_LSU.R /home/deakig/projects/metagenomics/data/fasta/ITS/S92_R2/ "*.\\.58" "*.\\.lsu" /home/deakig/projects/metagenomics/data/fasta/ITS/S92_R2.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_58Se_LSU.R /home/deakig/projects/metagenomics/data/fasta/ITS/S93_R2/ "*.\\.58" "*.\\.lsu" /home/deakig/projects/metagenomics/data/fasta/ITS/S93_R2.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_58Se_LSU.R /home/deakig/projects/metagenomics/data/fasta/ITS/S94_R2/ "*.\\.58" "*.\\.lsu" /home/deakig/projects/metagenomics/data/fasta/ITS/S94_R2.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_58Se_LSU.R /home/deakig/projects/metagenomics/data/fasta/ITS/S95_R2/ "*.\\.58" "*.\\.lsu" /home/deakig/projects/metagenomics/data/fasta/ITS/S95_R2.fa
	./ITS.sh /home/deakig/projects/metagenomics/rm_58Se_LSU.R /home/deakig/projects/metagenomics/data/fasta/ITS/S96_R2/ "*.\\.58" "*.\\.lsu" /home/deakig/projects/metagenomics/data/fasta/ITS/S96_R2.fa

### Remove chimeras
Using UNITE v 7.0 ITS database for chimeras (UCHIME reference dataset) https://unite.ut.ee/repository.php#uchime

	counter=91
	counter2=1
	for f in /home/deakig/projects/metagenomics/data/fasta/ITS/*.fa
	do 
	  if (( $counter2==1 ))
	  then
	    ./chimeras.sh $f /home/deakig/projects/metagenomics/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta S${counter}.${counter2}.cfree.fa /home/deakig/projects/metagenomics/data/fasta/de_chimerad/
	    counter2=2
	  else
	    ./chimeras.sh $f /home/deakig/projects/metagenomics/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta S${counter}.${counter2}.cfree.fa /home/deakig/projects/metagenomics/data/fasta/de_chimerad/
	    counter2=1
	    counter=$((counter+1));	
	  fi
	done

### Return merged common ITS1 and ITS2, unique ITS1 and unique ITS2

	./catfiles.pl S91.1.cfree.fa S91.2.cfree.fa S91
	./catfiles.pl S92.1.cfree.fa S92.2.cfree.fa S92
	./catfiles.pl S93.1.cfree.fa S93.2.cfree.fa S93
	./catfiles.pl S94.1.cfree.fa S94.2.cfree.fa S94
	./catfiles.pl S95.1.cfree.fa S95.2.cfree.fa S95
	./catfiles.pl S96.1.cfree.fa S96.2.cfree.fa S96

### OTU Picking and descriptive statistics
Multiple analyses were perfomed on:

1. Common and unique
2. Common ITS1 and ITS2
3. Unique ITS1
4. Unique ITS2

##### Common and unique (ITS1 and ITS2)
	cat *.fa > ITS.all.fa
	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/ITS/ITS.all.fa /home/deakig/projects/metagenomics/analysis/ITS_all_otus /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
	biom summarize-table -i ../analysis/ITS_all_otus/otu_table_mc2_w_tax.biom
	X=`biom summarize-table -i ../analysis/ITS_all_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
	./core_diversity.sh /home/deakig/projects/metagenomics/analysis/ITS_all_otus/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/analysis/ITS_all_cdout/ /home/deakig/projects/metagenomics/data/map.tsv . $X

##### Common ITS 
	cat S91.fa S92.fa S93.fa S94.fa S95.fa S96.fa > ITS.common.fa
	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/ITS/ITS.common.fa /home/deakig/projects/metagenomics/analysis/ITS_common_otus /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
	biom summarize-table -i ../analysis/ITS_common_otus/otu_table_mc2_w_tax.biom
	X=`biom summarize-table -i ../analysis/ITS_common_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
	./core_diversity.sh /home/deakig/projects/metagenomics/analysis/ITS_common_otus/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/analysis/ITS_common_cdout/ /home/deakig/projects/metagenomics/data/map.tsv . $X

##### Unique ITS1 only
	cat *r1* >ITS1.only.fa
	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/ITS/ITS1.only.fa /home/deakig/projects/metagenomics/analysis/ITS1_only_otus /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
	biom summarize-table -i ../analysis/ITS1_only_otus/otu_table_mc2_w_tax.biom
	X=`biom summarize-table -i ../analysis/ITS1_only_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
	./core_diversity.sh /home/deakig/projects/metagenomics/analysis/ITS1_only_otus/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/analysis/ITS1_only_cdout/ /home/deakig/projects/metagenomics/data/map.tsv . $X

##### Unique ITS2 only
	cat *r2* >ITS2.only.fa
	./pick_OTU.sh  /home/deakig/projects/metagenomics/data/fasta/ITS/ITS2.only.fa /home/deakig/projects/metagenomics/analysis/ITS2_only_otus /home/deakig/projects/metagenomics/scripts/params.txt /home/deakig/projects/metagenomics/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
	biom summarize-table -i ../analysis/ITS2_only_otus/otu_table_mc2_w_tax.biom
	X=`biom summarize-table -i ../analysis/ITS2_only_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
	./core_diversity.sh /home/deakig/projects/metagenomics/analysis/ITS2_only_otus/otu_table_mc2_w_tax.biom /home/deakig/projects/metagenomics/analysis/ITS2_only_cdout/ /home/deakig/projects/metagenomics/data/map.tsv . $X

### Statistical analysis
Requires a colData file describing condition (e.g. infected or uninfected) for each sample

As there were multiple OTU picking steps, multiple statistical analyses are necessary.

analysis.R biom_table "no. samples" median/geomean outfile
	
	Rscript analysis.R analysis/ITS_all_otus/otu_table_mc2_w_tax.biom colData median ITS.all.median.csv
	Rscript analysis.R analysis/ITS_common_otus/otu_table_mc2_w_tax.biom colData median ITS.median.csv
	Rscript analysis.R analysis/ITS1_only_otus/otu_table_mc2_w_tax.biom colData median ITS1.median.csv
	Rscript analysis.R analysis/ITS2_only_otus/otu_table_mc2_w_tax.biom colData median ITS2.median.csv
