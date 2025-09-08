#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J markadapt_rg
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

###################################### Parse the command inputs from wrapper script #########################

SAMPLE_ID=$1
FILENAME=$2
REFERENCE=$3
SAMPLE_TABLE=$4
READS_DIR=$5
OUTPUT_DIR=$6
INTERMEDIATES=$7
TASK_NUMBER=$8
READCOUNT=$9
CLUSTER=${10}
THREADS=${11}
ACCOUNT=${12}

####################################### load modules  ########################################################

# load modules
module load bioinfo-tools samtools/1.14 picard/2.23.4

###################### Get the reads and specify readgroups and some paths ###################################

# paths to reads from reads dir and sampletable
READS_1=${READS_DIR}/$(cat ${SAMPLE_TABLE} | awk -v sample_id=${SAMPLE_ID} ' $2 == sample_id ' | sed -n ${TASK_NUMBER}p | cut -f 4)
READS_2=${READS_DIR}/$(cat ${SAMPLE_TABLE} | awk -v sample_id=${SAMPLE_ID} ' $2 == sample_id ' | sed -n ${TASK_NUMBER}p | cut -f 5)

# make a "run name" that discriminates between the different runs (if split sample), this is what we'll use for storing tempfiles
RUN_NAME=${FILENAME}_$(basename ${READS_1%%.fastq*})

# add read group data
#read group ID
ID=$(cat ${SAMPLE_TABLE} | grep $(basename ${READS_1}) | cut -f 1  | cut -d "_" -f 1-2)
#platform unit
PU=$(cat ${SAMPLE_TABLE} | grep $(basename ${READS_1}) | cut -f 1)
#sample name
SM=${FILENAME}
#sequencing platform
PL=ILLUMINA
#library name
LB=${FILENAME}
#insert size
PI=350

# in the next step I'll map this data, sending of the mapping script to the queue already here, with dependency to wait for this one to finish
sbatch --dependency=afterok:${SLURM_JOB_ID} -M ${CLUSTER} -J ${SAMPLE_ID}.map.${TASK_NUMBER} -A ${ACCOUNT} 2_mapping.sh \
	${INTERMEDIATES}/${RUN_NAME}.fastq.gz \
	${INTERMEDIATES}/${RUN_NAME}.u.bam \
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

# Let the log now what we're doing
echo $(date)": Task "${TASK_NUMBER}"of"${READCOUNT} " started processing sample "${FILENAME} > ${OUTPUT_DIR}/${FILENAME}.log

# check if we already have some of the starting files, in which case we'll skip those steps

if [ ! -f ${INTERMEDIATES}/${RUN_NAME}.u.bam ]; then
echo "adding readgroup data and marking adapters to readpair:" >> ${OUTPUT_DIR}/${FILENAME}.log
echo ${READS_1} >> ${OUTPUT_DIR}/${FILENAME}.log
echo ${READS_2} >> ${OUTPUT_DIR}/${FILENAME}.log

################################## Run the software ###########################################
if [ ! -f ${TASK_NUMBER} -eq 1 ]; then
echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "Modules/versions used in the preprocessing reads script: \n" >> ${OUTPUT_DIR}/${FILENAME}.log
module list >> ${OUTPUT_DIR}/${FILENAME}.log

echo -e "FastqToSam command run as follows:\njava -jar ${PICARD_HOME}/picard.jar FastqToSam \\
\tF1=${READS_1} \\
\tF2=${READS_2} \\
\tSM=${SM} \\
\tRG=${ID} \\
\tPU=${PU} \\
\tPL=${PL} \\
\tLB=${LB} \\
\tPI=${PI} \\
\tO=${INTERMEDIATES}/${RUN_NAME}.u.bam" >> ${OUTPUT_DIR}/${FILENAME}.log
fi

java -jar ${PICARD_HOME}/picard.jar FastqToSam \
	F1=${READS_1} \
	F2=${READS_2} \
	SM=${SM} \
	RG=${ID} \
	PU=${PU} \
	PL=${PL} \
	LB=${LB} \
	PI=${PI} \
	O=${INTERMEDIATES}/${RUN_NAME}.u.bam
else
echo "Unmapped bam already present, skipping this step"
fi
if [ ! -f ${INTERMEDIATES}/mi.${RUN_NAME}.u.bam ]; then
if [ ! -f ${INTERMEDIATES}/${RUN_NAME}.fastq.gz ]; then
#next we use MarkIlluminaAdapters to mark adapter content in the reads, creating one intermediate bam-file and one textfile containing info on adapter content
if [ ${TASK_NUMBER} -eq 1 ]; then
echo -e "MarkIlluminaAdapters command run as:\n
java -jar ${PICARD_HOME}/picard.jar MarkIlluminaAdapters \\
\tI=${INTERMEDIATES}/${RUN_NAME}.u.bam \\
\tO=${INTERMEDIATES}/mi.${RUN_NAME}.u.bam \\
\tM=${OUTPUT_DIR}/mi.${RUN_NAME}.metrix.txt\n" >> ${OUTPUT_DIR}/${FILENAME}.log
fi
java -jar ${PICARD_HOME}/picard.jar MarkIlluminaAdapters \
	I=${INTERMEDIATES}/${RUN_NAME}.u.bam \
	O=${INTERMEDIATES}/mi.${RUN_NAME}.u.bam \
	M=${OUTPUT_DIR}/mi.${RUN_NAME}.metrix.txt

else
echo "Marked fastq already present, skipping this step" >> ${OUTPUT_DIR}/${FILENAME}.log
fi
else
echo "Marked bam already present, skipping this step" >> ${OUTPUT_DIR}/${FILENAME}.log
fi
echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log

if [ ! -f ${INTERMEDIATES}/${RUN_NAME}.fastq.gz ]; then
#then revert this file back to an interleaved fastq file, which we'll again save to the tempdir and use for mapping in the next step
if [ ${TASK_NUMBER} -eq 1 ]; then
echo -e "\nSamToFastq command run as: \n
java -jar ${PICARD_HOME}/picard.jar SamToFastq \\
\tI=${INTERMEDIATES}/mi.${RUN_NAME}.u.bam \\
\tFASTQ=${INTERMEDIATES}/${RUN_NAME}.fastq.gz \\
\tCLIPPING_ATTRIBUTE=XT \\
\tCLIPPING_ACTION=2 \\
\tINTERLEAVE=true \\
\tNON_PF=true" >> ${OUTPUT_DIR}/${FILENAME}.log
fi
java -jar ${PICARD_HOME}/picard.jar SamToFastq \
	I=${INTERMEDIATES}/mi.${RUN_NAME}.u.bam \
	FASTQ=${INTERMEDIATES}/${RUN_NAME}.fastq.gz \
	CLIPPING_ATTRIBUTE=XT \
	CLIPPING_ACTION=2 \
	INTERLEAVE=true \
	NON_PF=true 
# remove the mi-bam
rm ${INTERMEDIATES}/mi.${RUN_NAME}.u.bam

else
echo "Transformed fastq file already present, skipping this step."
fi
echo $(date)": adapters marked and readgroup data added in task " ${TASK_NUMBER} " of " ${READCOUNT} ", now proceeding to mapping." >> ${OUTPUT_DIR}/${FILENAME}.log

if [ ${TASK_NUMBER} -eq 1 ]; then
echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
fi

