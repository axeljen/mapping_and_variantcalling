#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5-00:00:00
#SBATCH -J genotype
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/genotype-%x-%j.out
#SBATCH -e ./logs/genotype-%x-%j.error

#load modules
module load bioinfo-tools GATK/4.3.0.0

# parse command line inputs
REFERENCE=${1}
ARGFILE=${2}
OUTPUT_DIR=${3}
CHROM=${4}
START=${5}
END=${6}
RUN_NAME=${7}

mkdir -p ${OUTPUT_DIR}

#then put together the output string
OUT_VCF=${OUTPUT_DIR}/${RUN_NAME}_${CHROM}_${START}_${END}.vcf.gz

#first combine all the gvfs on temp-storage
gatk --java-options "-Xmx6G" CombineGVCFs \
	-R ${REFERENCE} \
	-O ${SNIC_TMP}/${CHROM}.g.vcf \
	-L ${CHROM}:${START}-${END} \
	--arguments_file ${ARGFILE}

#then genotype
gatk --java-options "-Xmx6G" GenotypeGVCFs \
	-R ${REFERENCE} \
	-V ${SNIC_TMP}/${CHROM}.g.vcf \
	-L ${CHROM}:${START}-${END} \
	-all-sites \
	-O ${OUT_VCF}
