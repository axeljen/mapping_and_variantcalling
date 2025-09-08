#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH -J mapping
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

############################# load modules ##############################

module load bioinfo-tools samtools/1.14 picard/2.23.4 bwa/0.7.17

################################ Parse command inputs from 1_preprocess_reads.sh #########################
	
READS=${1}
UBAM=${2}
SAMPLE_ID=${3}
FILENAME=${4}
REFERENCE=${5}
SAMPLE_TABLE=${6}
READS_DIR=${7}
OUTPUT_DIR=${8}
INTERMEDIATES=${9}
TASK_NUMBER=${10}
READCOUNT=${11}
RUN_NAME=${12}
CLUSTER=${13}
THREADS=${14}
ACCOUNT=${15}


# send off the merging script to queue with dependency to wait for this one
sbatch --dependency=afterok:${SLURM_JOB_ID} -M ${CLUSTER} -J ${SAMPLE_ID}.mergebams.${TASK_NUMBER} -A ${ACCOUNT} 3_mergebams.sh \
	${INTERMEDIATES}/${RUN_NAME}.u.bam \
	${INTERMEDIATES}/sorted.${RUN_NAME}.bam \
	${SAMPLE_ID} \
	${FILENAME} \
	${REFERENCE} \
	${SAMPLE_TABLE} \
	${READS_DIR} \
	${OUTPUT_DIR} \
	${INTERMEDIATES} \
	${TASK_NUMBER} \
	${READCOUNT} \
	${RUN_NAME} \
	${CLUSTER} \
	${THREADS} \
	${ACCOUNT}

################################## Map ######################################

echo $(date) ": This is task "${TASK_NUMBER}"/"${READCOUNT}" for sample "${SAMPLE_ID}", now mapping..." >> ${OUTPUT_DIR}/${FILENAME}.log

if [ ! -f ${INTERMEDIATES}/sorted.${RUN_NAME}.bam ]; then
# mapping the interleaved, processed reads from the processing script, sorting this and putting it in a temporary file

# First, if this is task 1, log to the "master log"
if [ ${TASK_NUMBER} -eq 1 ]; then
echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "Modules/versions used in the mapping/preprocessing script: \n" >> ${OUTPUT_DIR}/${FILENAME}.log
module list >> ${OUTPUT_DIR}/${FILENAME}.log

echo -e "\nMapping command run as:\n
bwa mem -M -t ${THREADS} -p ${REFERENCE} interleaved.reads.fq | \\
java -Xmx120G -jar ${PICARD_HOME}/picard.jar SortSam \\
\tINPUT=/dev/stdin \\
\tOUTPUT=intermediates/sorted.filename.bam \\
\tSORT_ORDER=coordinate \\
\tCREATE_INDEX=true" >> ${OUTPUT_DIR}/${FILENAME}.log
fi
bwa mem -M -t ${THREADS} -p ${REFERENCE} ${READS} | \
java -Xmx120G -jar ${PICARD_HOME}/picard.jar SortSam \
	INPUT=/dev/stdin \
	OUTPUT=${INTERMEDIATES}/sorted.${RUN_NAME}.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true

echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log

else
echo "Mapped bam already present, skipping this step" >> ${OUTPUT_DIR}/${FILENAME}.log
fi
# last check if mapping was successful, in which case we can delete the input fastqreads
outbam=${INTERMEDIATES}/sorted.${RUN_NAME}.bam

samtools quickcheck ${outbam} && rm ${READS} || echo "Something went wrong with the mapping it seems, not deleting any intermediates for now." >> ${OUTPUT_DIR}/${FILENAME}.log

