#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J mergebams
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

# this job runs two tools: mergebamalignment from picard

#load modules
module load bioinfo-tools samtools/1.14 picard/2.23.4

# parse all the command inputs from previous job
UBAM=${1}
MAPPED_BAM=${2}
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

echo $(date)": this is MergeBamAlignment for task "${TASK_NUMBER} "of" ${READCOUNT}", retaining readgroup info from unmapped bam, merging it with the merged alignment." >> ${OUTPUT_DIR}/${FILENAME}.log


#Now we use MergeBamAlignment to retain the metadata and add some additional info to our mapped bamfile,
#using the unmapped bam we created in step one

# check if this step already has been performed
samtools quickcheck ${INTERMEDIATES}/mba.${RUN_NAME}.bam && MBA_EXISTS=TRUE

if [ -z ${MBA_EXISTS} ]; then
if [ ${TASK_NUMBER} -eq 1 ]; then
echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "Modules/versions used in the MergeBamAlignment script: \n" >> ${OUTPUT_DIR}/${FILENAME}.log
module list >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "MergeBamAlignment command run as: \n
java -Xmx4G -jar ${PICARD_HOME}/picard.jar MergeBamAlignment \
\tR=${REFERENCE} \
\tUNMAPPED_BAM=${UBAM} \
\tALIGNED_BAM=${MAPPED_BAM} \
\tO=${INTERMEDIATES}/mba.${RUN_NAME}.bam \
\tCREATE_INDEX=true \
\tADD_MATE_CIGAR=true \
\tCLIP_ADAPTERS=false \
\tCLIP_OVERLAPPING_READS=true \
\tINCLUDE_SECONDARY_ALIGNMENTS=true \
\tMAX_INSERTIONS_OR_DELETIONS=-1 \
\tPRIMARY_ALIGNMENT_STRATEGY=MostDistant \
\tATTRIBUTES_TO_RETAIN=XS" >> ${OUTPUT_DIR}/${FILENAME}.log

fi
java -Xmx4G -jar ${PICARD_HOME}/picard.jar MergeBamAlignment \
	R=${REFERENCE} \
	UNMAPPED_BAM=${UBAM} \
	ALIGNED_BAM=${MAPPED_BAM} \
	O=${INTERMEDIATES}/mba.${RUN_NAME}.bam \
	CREATE_INDEX=true \
	ADD_MATE_CIGAR=true \
	CLIP_ADAPTERS=false \
	CLIP_OVERLAPPING_READS=true \
	INCLUDE_SECONDARY_ALIGNMENTS=true \
	MAX_INSERTIONS_OR_DELETIONS=-1 \
	PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
	ATTRIBUTES_TO_RETAIN=XS

echo ${RUN_NAME} >> ${INTERMEDIATES}/merged_bams.txt

else
echo " mba for this run alread performed, skipping ahead.." >> ${OUTPUT_DIR}/${FILENAME}.log
fi
# if this step was succesful, we can delete the input bams
outbam=${INTERMEDIATES}/mba.${RUN_NAME}.bam

samtools quickcheck ${outbam} && rm ${UBAM} ${MAPPED_BAM} || echo "Error: something went wrong in the merging step, not deleting intermediates right now." >> ${OUTPUT_DIR}/${FILENAME}.log

# next we need to check if this sample was sequenced in multiple lanes - if so these should be merged at this step
# if n is more than 1, we take a small detour to check if this sample is ready for merging
if [ ${READCOUNT} -gt 1 ];
	then
	# ok so now we need to check if this is the last bam of all for this sample to finish, and if so merge these, otherwise skip this step
	count=$(cat ${INTERMEDIATES}/merged_bams.txt | wc -l)
	#echo "There are " ${indeces} " bamfiles in the directory." >> ${OUTDIR}/${FILENAME}.log
	# count and see if we should proceed
	if [ ${count} -eq ${READCOUNT} ];
		then
		echo "All" ${READCOUNT} " bams seems to be processed for this sample. Continuing with merging these." >> ${OUTPUT_DIR}/${FILENAME}.log
		# put all the bamfiles in a list
		find ${INTERMEDIATES} -name "*mba*.bam" > ${INTERMEDIATES}/bamfiles.txt
		# send it off to samtools merge
		sbatch -M ${CLUSTER} -J ${SAMPLE_ID}.merge -A ${ACCOUNT} samtools_merge.sh \
			${INTERMEDIATES}/bamfiles.txt \
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
	fi
else
	# if n = 1, we just proceed to markdups directly
	echo "This sample was sequenced on a single lane. Sending it off to MarkDuplicates." >> ${OUTPUT_DIR}/${FILENAME}.log
	sbatch -M ${CLUSTER} -J ${SAMPLE_ID}.markdups -A ${ACCOUNT} 5_markdups.sh \
		${INTERMEDIATES}/mba.${RUN_NAME}.bam \
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
fi
echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log

