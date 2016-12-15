#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

for thing in "$@";do
	echo removing $thing
	rm -r $thing
done