# Mapping

This directory contains the scripts required to go from raw fastq reads to mapped and deduplicated bam files, ready for variant calling.

### Dependencies

The scripts in this directory uses the following software:
- samtools/1.14
- bwa/0.7.17
- picard/2.23.4

These are available as preinstalled modules on the cluster for which this pipeline was developed (Uppmax).

### Running the pipeline

1) Set up the sample table, should contain one row per sample/lane, containing the following five tab-separated columns:
Sequencing_ID	Sample_id	Taxon	forward_reads	reverse_reads

Only the basename of the reads should be provided, as the path to their directory is specified in the start_mapping.sh script. The sequencing ID is really only relevant if you have multiple lanes per sample, otherwise it can be the same as the sample_id. Taxon column is used for naming output files.

2) Make sure that the reference genome is indexed with bwa, samtools and picard:
```
bwa index reference.fasta
samtools faidx reference.fasta
picard CreateSequenceDictionary R=reference.fasta
```

3) Modify the parameters and paths in the first section of the start_mapping.sh script (nothing below line 32 should need to be changed).

4) The start_mapping.sh takes the sample_id as command line argument, so to start the mapping for all samples in the sample table, you can use a loop like this:
```
for sample in $(cut -f2 sample_table_example.txt | tail -n +2); do
    bash start_mapping.sh $sample
done
```
