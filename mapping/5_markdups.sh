#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-20:00:00
#SBATCH -J dedup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

# load modules
module load bioinfo-tools samtools/1.14 picard/2.23.4

##################### Parse inputs from previous step
BAMFILE_IN=${1}
SAMPLE_ID=${2}
FILENAME=${3}
REFERENCE=${4}
SAMPLE_TABLE=${5}
READS_DIR=${6}
OUTPUT_DIR=${7}
INTERMEDIATES=${8}
READCOUNT=${9}
CLUSTER=${10}
THREADS=${11}
ACCOUNT=${12}

# using string alterations to save the dedupped bamfile in the same directory, prefixed md
BAMFILE_OUT=${OUTPUT_DIR}/md.${FILENAME}.bam

# we'll also save the metrix file produced by markduplicates there, replacing the ".bam" suffix to "_metrix.txt"
METRIX=${BAMFILE_OUT/.bam/_metrix.txt}

echo -e "\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
echo $(date)": This is MarkDuplicates, processing this bamfile: " >> ${OUTPUT_DIR}/${FILENAME}.log
echo ${BAMFILE_IN} >> ${OUTPUT_DIR}/${FILENAME}.log
echo "Final output file will be stored here: " >> ${OUTPUT_DIR}/${FILENAME}.log
echo ${OUTDIR}/md.${FILENAME}.bam >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "Modules/versions used in the MarkDuplicates script: \n" >> ${OUTPUT_DIR}/${FILENAME}.log
module list >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "MarkDuplicates command run as:\n" >> ${OUTPUT_DIR}/${FILENAME}.log
echo -e "java -Xmx5g -jar ${PICARD_HOME}/picard.jar MarkDuplicates \\
        I=${BAMFILE_IN} \\
        O=${BAMFILE_OUT} \\
        M=${METRIX} \\
        ASSUME_SORT_ORDER=coordinate \\
	MAX_RECORDS_IN_RAM=750000" >> ${OUTPUT_DIR}/${FILENAME}.log

# submit the qualimap to run afterok
sbatch -M ${CLUSTER} -A ${ACCOUNT} --dependency=afterok:${SLURM_JOB_ID} 6_qualimap.sh ${BAMFILE_OUT} ${REFERENCE} ${CLUSTER}

# then run the mark duplicates tool
java -Xmx5g -jar ${PICARD_HOME}/picard.jar MarkDuplicates \
        I=${BAMFILE_IN} \
        O=${BAMFILE_OUT} \
        M=${METRIX} \
        ASSUME_SORT_ORDER=coordinate \
	MAX_RECORDS_IN_RAM=750000

# index this one
samtools index ${BAMFILE_OUT}

# check if succesful
outbam=${OUTDIR}/md.${FILENAME}.bam

samtools quickcheck ${BAMFILE_OUT} && echo "All done. Final bamfile, ready for variant calling, is here:"; \
        echo ${BAMFILE_OUT}  >> ${OUTPUT_DIR}/${FILENAME}.log

samtools quickcheck ${BAMFILE_OUT} && ALL_OK=TRUE

if [ ! -z ${ALL_OK} ]; then
echo "Removing intermediate directory." >> ${OUTPUT_DIR}/${FILENAME}.log
rm -r ${INTERMEDIATES}
fi

echo -e "\nAll done.\n##########################################################\n" >> ${OUTPUT_DIR}/${FILENAME}.log
