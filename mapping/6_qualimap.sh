#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-20:00:00
#SBATCH -J qualimap
#SBATCH --mail-type=fail
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/qualimap-%j.out
#SBATCH -e ./logs/qualimap-%j.error

#load modules
module load bioinfo-tools QualiMap

#some bug in qualimap has forced me to add this
unset DISPLAY

#input bam as command line input
BAMFILE_IN=$1

# directory to store output
OUTDIR=${BAMFILE_IN%}.mapstats
mkdir -p ${OUTDIR}

#then run it, it's quite memory heavy, so we give it a full node and java memory of 140G
qualimap bamqc --java-mem-size=120G \
	-bam ${BAMFILE_IN} \
	-nt 20 \
	-outformat PDF:HTML \
	-outdir ${OUTDIR}
