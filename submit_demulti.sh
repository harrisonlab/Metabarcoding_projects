#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

for SCRIPT_DIR; do true; done

${SCRIPT_DIR}/demulti_v2.pl $@


