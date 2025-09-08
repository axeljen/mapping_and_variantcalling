######################################### Wrapper script to run haplotypecaller single sample, parellelized by splitting up reference genome

######################################## Modify these paths before executing
# Reference genome
REFERENCE=../testfiles/reference/GCF_003339765.1_Mmul_10_NC_041770.1.fna

# path to a file that lists chromosomes/scaffolds and their lengths, variants will only be called on these regions.
CHROMOSOMES=reference_chromosomes.txt

# Name of this run (will be used for naming the parent directory and final vcf)
RUN_NAME=MMUL_10

# Where to store the output gvcf, will be stored in a sample specific directory within this one, and will use this dir for some bedfiles etc as well
OUTPUT_DIR=../variants

# For parallelization, split the genome into chunks of this size (in bp)
INTERVAL_SIZE=20000000

# Directory where to our bamfiles are stored. All and only the bamfiles that we wanna call variants from should be in this or any child directory, as the pipelien is set up.
BAMDIR=../mappings/

# account to use on the cluster
ACCOUNT=uppmax2025-2-295

# cluster to run calling in
CLUSTER=rackham

##################################################################################
mkdir -p ${OUTPUT_DIR}


# use the reference_intervals.py script to create a file with intervals of this size (will create reference_intervals.txt)
python3 reference_intervals.py ${CHROMOSOMES} ${INTERVAL_SIZE} 

INTERVALS=reference_intervals.txt

# make a list of bamfiles for which to call variants
find ${BAMDIR} -name "*.bam" > bamfiles.txt
BAMFILES=bamfiles.txt

# last we will submit an array job for each interval, where each array should have as many tasks as there are samples (cramfiles)
while read interval
do
CHROM=$(echo ${interval} | cut -d " " -f 1)
START=$(echo ${interval} | cut -d " " -f 2)
END=$(echo ${interval} | cut -d " " -f 3)
NTASKS=$(cat ${BAMFILES} | wc -l)
sbatch -M ${CLUSTER} -A ${ACCOUNT} --array=1-${NTASKS} -J ${RUN_NAME}.hapcaller 1_hapcaller.sh \
	${REFERENCE} \
	${BAMFILES} \
	${CHROM} \
	${START} \
	${END} \
	${OUTPUT_DIR} \
	${RUN_NAME} \
    ${ACCOUNT}
done < ${INTERVALS}
