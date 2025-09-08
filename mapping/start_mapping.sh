

# Wrapper/started for mapping and processing pipeline following gatk workflow

# Axel Jensen, 2022

######################################################################################################################################
# a couple of parameters that should be modified here:

# specify what cluster to run on
CLUSTER=rackham

# number of maximum threads to use in scripts that multithread
THREADS=20

# project id for submitting jobs
ACCOUNT=uppmax2025-2-295

# will need a sample table that we will use for getting the appropriate reads, naming the files etc. should have five, tab-separated, colummns:sequencing_id,sample_id,taxon,forward_reads(just basename is ok),reverse_reads
SAMPLE_TABLE=sample_table_example.txt

# all read files for the project should be stored (or linked to) in a single directory, specify this one here
READS_DIR=../testfiles/reads

# where to store the final output file
OUTPUT_DIR=../mappings

# make this dir if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Reference genome for aligning the reads to (should be indexed and ready)
REFERENCE=../testfiles/reference/GCF_003339765.1_Mmul_10_NC_041770.1.fna

#### the following scripts needs to be located in the same directory as this wrapper, and will be run in the following order
# 1_preprocess_reads.sh
# 2_mapping.sh
# 3_mergebams.sh
# 4_samtools_merge.sh
# 5_markdups.sh
# 6_qualimap.sh

#####################################################################################################################################

# take the sample id as command input
SAMPLE_ID=$1

set -e

# from this id, fetch all the forward reads and count how many read pairs it has
READCOUNT=$(cat ${SAMPLE_TABLE} | awk -v sample_id=${SAMPLE_ID} ' $2 == sample_id ' | cut -f 4 | wc -l)

# get the sample name/filename, constructed from sample id + taxon
TAXON=$(cat ${SAMPLE_TABLE} | awk -v sample_id=${SAMPLE_ID} ' $2 == sample_id ' | head -n 1 | cut -f 3)
# Construct the filename using sample ID and taxon
FILENAME=${SAMPLE_ID}_${TAXON}

# a directory that we'll use for storing some intermediates along the way, and then delete at the end
INTERMEDIATES=intermediates.${FILENAME}
mkdir -p ${INTERMEDIATES}

# start a little custom logfile where we can output some progression as the pipeline goes ahead
LOGFILE=${OUTPUT_DIR}/${FILENAME}.log

echo $(date)": Starting the processing and mapping pipeline for sample " ${FILENAME}"." > ${LOGFILE}
if [ ${READCOUNT} -gt 1 ];
then
echo -e "\nThis sample is split across "${READCOUNT}" read pairs.\n" >> ${LOGFILE}
else
echo -e "\nThis sample contains a single read pair.\n" >> ${LOGFILE}
fi
echo -e "Final output will be stored here: \n" >> ${LOGFILE}
echo -e ${OUTPUT_DIR}/${FILENAME}.bam "\n" >> ${LOGFILE}

echo -e "\n########################################################################\n"  >> ${LOGFILE}

######################################## Print all variables for checking
echo -e "Variables and paths set as follows: \\
CLUSTER=${CLUSTER} \\
THREADS=${THREADS} \\
INTERMEDIATES=${INTERMEDIATES} \\
SAMPLE_ID=${SAMPLE_ID} \\
FILENAME=${FILENAME} \\
REFERENCE=${REFERENCE} \\
SAMPLE_TABLE=${SAMPLE_TABLE} \\
READS_DIR=${READS_DIR} \\
OUTPUT_DIR=${OUTPUT_DIR} \\
READCOUNT=${READCOUNT}" >> ${LOGFILE}

echo -e "\n########################################################################\n"  >> ${LOGFILE}

# Check at what stage to start this job, this is a quick and dirty check, if the files exist but the last script failed, it will try to start from there anyway. So a bit of manual checking for failing jobs might be needed.

# if the "final" mba file is in the INTERMEDIATES, start directly at dedupping
if [ -f ${INTERMEDIATES}/mba.${FILENAME}.bam ]; then
echo -e "Found a merged bamfile here, staring with MarkDuplicates. \n" >> ${LOGFILE}.log
samtools quickcheck ${INTERMEDIATES}/mba.${FILENAME}.bam && \
sbatch -M ${CLUSTER} -J ${SAMPLE_ID}.dedup -A ${ACCOUNT} 5_markdups.sh \
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

# if there are any mba files in the INTERMEDIATES, we start at the merging step
elif [ $(find ${INTERMEDIATES} -name "*mba*.bam" | wc -l) -gt 0 ]; then
echo -e "Found merged bamfiles in this directory, starting at this step already. \n" >> ${LOGFILE}.log
# if READCOUNT is more than 1, we take a small detour to check if this sample is ready for merging
if [ ${READCOUNT} -gt 1 ];
	then
	# ok so now we need to check if this is the last bam of all for this sample to finish, and if so merge these, otherwise skip this step
	indeces=$(find ${INTERMEDIATES} -name "*mba.*.bai" | wc -l)
	#echo "There are " ${indeces} " bamfiles in the directory." >> ${OUTDIR}/${FILENAME}.log
	# count and see if we should proceed
	if [ ${indeces} -eq ${READCOUNT} ];
		then
		echo -e "All" ${READCOUNT} " bams seems to be processed for this sample. Continuing with merging these. \n" >> ${LOGFILE}
		# put all the bamfiles in a list
		find ${INTERMEDIATES} -name "*mba*.bam" > ${INTERMEDIATES}/bamfiles.txt
		# send it off to samtools merge
		sbatch -M ${CLUSTER} -J ${SAMPLE_ID}.merge -A ${ACCOUNT} 4_samtools_merge.sh \
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
	# if READCOUNT = 1, we just proceed to markdups directly
	if [ ${READCOUNT} -eq 1 ];
	then
	echo -e "This sample was sequenced on a single lane. Sending it off to MarkDuplicates.\n" >> ${LOGFILE}
	# get the mba file
	BAMFILE=$(find ${INTERMEDIATES} -name "mba.*.bam")
	samtools quickcheck ${BAMFILE} && \
	sbatch -M ${CLUSTER} -A ${ACCOUNT} 5_markdups.sh \
		${BAMFILE} \
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
fi

# else check if we have mapped bams, in which case we don't need to map and start on mergebams
elif [ $(find ${INTERMEDIATES} -name "sorted.*.bam" | wc -l) -gt 0 ]; then
echo -e "Found mapped bamfiles in this directory, starting with merging the bamfiles for each run. \n" >> ${LOGFILE}
for i in $(seq 1 $(find ${INTERMEDIATES} -name "sorted.*.bam" | wc -l));
do
bam=$(find ${INTERMEDIATES} -name "sorted.*.bam" | sed -n ${i}p)
TASK_NUMBER=${i}
RUN_NAME=$(basename ${bam} | cut -d "." -f 2 | rev | cut -d "." -f 2- | rev)
samtools quickcheck ${bam} && sbatch -M ${CLUSTER} -A ${ACCOUNT} -J ${SAMPLE_ID}.mergebams.${TASK_NUMBER} 3_mergebams.sh \
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
done
# else check if preprocessing of reads finished and in that case start the mapping
elif [ $(find ${INTERMEDIATES} -name "*.fastq.gz" | wc -l) -gt 0 ]; then
echo -e "Found processed readfiles in this directory already, starting at the mapping stage. \n" >> ${LOGFILE}
for i in $(seq 1 $(find ${INTERMEDIATES} -name "*.fastq.gz" | wc -l));
do
READS=$(find ${INTERMEDIATES} -name "*.fastq.gz" | sed -n ${i}p)
TASK_NUMBER=${i}
RUN_NAME=$(basename ${READS%%.fastq.gz})
sbatch -M ${CLUSTER} -J ${SAMPLE_ID}.map.${TASK_NUMBER} -A ${ACCOUNT} 2_mapping.sh \
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
done
# if nothing of the above, start the read processing
else
# next we submit the job for processing and mapping, a sepearate job for each read pair
for i in $(seq 1 ${READCOUNT});
do
TASK_NUMBER=${i}
sbatch -M ${CLUSTER} -J ${SAMPLE_ID}.mi.rg.${i} -A ${ACCOUNT} 1_preprocess_reads.sh \
	${SAMPLE_ID} \
	${FILENAME} \
	${REFERENCE} \
	${SAMPLE_TABLE} \
	${READS_DIR} \
	${OUTPUT_DIR} \
	${INTERMEDIATES} \
	${TASK_NUMBER} \
	${READCOUNT} \
	${CLUSTER} \
	${THREADS} \
	${ACCOUNT}
done
fi

echo -e "Sent of ${READCOUNT} readpairs for prerprocessing. \n" >> ${LOGFILE}