#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

R1=$1
ADAPTERS=$2

R2=$(echo $R1|sed 's/_R1_/_R2_/')

usearch8.1 -search_oligodb $R1 -db $ADAPTERS -strand both -userout ${R1}.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
usearch8.1 -search_oligodb $R2 -db $ADAPTERS -strand both -userout ${R2}.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
cat ${R1}.txt ${R2}.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' $R1
cat ${R1}.txt ${R2}.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' $R2
