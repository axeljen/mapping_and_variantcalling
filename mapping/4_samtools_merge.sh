#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 1-00:00:00
#SBATCH -J samtools_merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

# load modules
module load bioinfo-tools samtools/1.14


# parse inputs
BAMLIST=${1}
SAMPLE_ID=${2}
FILENAME=${3}
REFERENCE=${4}
SAMPLE_TABLE=${5}
READS_DIR=${6}
OUTPUT_DIR=${7}
INTERMEDIATES=${8}
#TASK_NUMBER=${9}
READCOUNT=${9}
CLUSTER=${10}
THREADS=${11}
ACCOUNT=${12}

echo ${CLUSTER}

echo "This is samtools merge merging the following bamfiles: " >> ${OUTPUT_DIR}/${FILENAME}.log
cat ${BAMLIST} >> ${OUTPUT_DIR}/${FILENAME}.log

echo "Storing the merged bamfile to " ${INTERMEDIATES}/mba.${FILENAME}.bam >> ${OUTPUT_DIR}/${FILENAME}.log

# submit markdups but make it wait for this one
sbatch --dependency=afterok:${SLURM_JOB_ID} -M ${CLUSTER} -J ${SAMPLE_ID}.dedup 5_markdups.sh \
	${INTERMEDIATES}/mba.${FILENAME}.bam \
	${SAMPLE_ID} \
	${FILENAME} \
	${REFERENCE} \
	${SAMPLE_TABLE} \
	${READS_DIR} \
	${OUTPUT_DIR} \
	${INTERMEDIATES} \
	${READCOUNT} \
	${CLUSTER} \
	${THREADS} \
	${ACCOUNT}

echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "Modules/versions used in the samtools merge script: \n" >> ${OUTPUT_DIR}/${FILENAME}.log
module list >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "\n samtools merge script run as:\n"
echo -e "samtools merge -@ 10 -b ${BAMLIST} ${INTERMEDIATES}/mba.${FILENAME}.bam\n"
# merge these bams
samtools merge -@ 10 -b ${BAMLIST} ${INTERMEDIATES}/mba.${FILENAME}.bam

# index it
samtools index ${INTERMEDIATES}/mba.${FILENAME}.bam

# check if succesful and in that case delete intermediates and proceed
outbam=${INTERMEDIATES}/mba.${FILENAME}.bam

samtools quickcheck ${outbam} && ALL_OK=TRUE

if [ ! -z ${ALL_OK} ];
then
echo "Merging succesful, deleting intermediate bams:" >> ${OUTPUT_DIR}/${FILENAME}.log
for bam in $(cat ${BAMLIST})
do

cat ${BAMLIST} >> ${OUTPUT_DIR}/${FILENAME}.log
# check inbam just as a precaution, not to delete anything we don't want to
samtools quickcheck ${bam} && rm ${bam}
done
fi

echo "Now sent of to MarkDuplicates." >> ${OUTPUT_DIR}/${FILENAME}.log

echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
