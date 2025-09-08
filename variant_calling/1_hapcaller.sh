#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH -J hapcaller
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/hapcaller-%x-%j.out
#SBATCH -e ./logs/hapcaller-%x-%j.error

# load modules
module load bioinfo-tools GATK/4.3.0.0 samtools/1.14

# specify what reference we've mapped to
REFERENCE=${1}
BAMFILES=${2}
CHROM=${3}
START=${4}
END=${5}
OUTPUT_DIR=${6}
RUN_NAME=${7}
ACCOUNT=${8}

# use the task number to get the bamfile
BAMFILE_IN=$(cat ${BAMFILES} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# if it's not indexed, fix
if [ ! -f ${BAMFILE_IN}.crai ]; then
samtools index ${BAMFILE_IN}
fi

# make a filename
FILENAME=$(basename ${BAMFILE_IN%%.bam})
# if it's a cram
FILENAME=${FILENAME%%.cram}
# add chrom, start and end
FILENAME=${FILENAME}_${CHROM}_${START}_${END}

# store gvcfs in a dedicated dir
GVCFs=${OUTPUT_DIR}/GVCFs

#first, if there is not yet a directory for the current chromosome in outdir, we create this
mkdir -p ${GVCFs}/${CHROM}

# make an outfile path
OUTFILE=${GVCFs}/${CHROM}/${FILENAME}.g.vcf.gz

# echo this to an argfile that we'll use for combining and genotyping
mkdir -p tmp.argfiles
ARGFILE=tmp.argfiles/${SLURM_ARRAY_JOB_ID}.args.txt

echo "-V" ${OUTFILE} >> ${ARGFILE}
        
# if this is the last task, submit the combine and genotype job
if [ ${SLURM_ARRAY_TASK_ID} -eq ${SLURM_ARRAY_TASK_MAX} ];
then
sbatch -M ${CLUSTER} -A ${ACCOUNT} -J ${RUN_NAME}.genotype --dependency=afterok:${SLURM_JOB_ID} 2_combine_genotype.sh \
    ${REFERENCE} \
    ${ARGFILE} \
    ${OUTPUT_DIR} \
    ${CHROM} \
    ${START} \
    ${END} \
    ${RUN_NAME}
fi

# run hapcaller here
gatk --java-options "-Xmx6G" HaplotypeCaller \
        -I ${BAMFILE_IN} \
        -O ${OUTFILE} \
        -R ${REFERENCE} \
		-L ${CHROM}:${START}-${END} \
        -ERC GVCF
