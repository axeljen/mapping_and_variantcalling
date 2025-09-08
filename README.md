# Short read variant calling pipeline

This is a semi-automated pipeline to process and map short read sequencing data to a reference genome, and call variants. The final output is a set of jointly genotyped VCF files.

The pipeline was written for a specific HPC cluster, and will need modifications to run on other systems.

The pipeline in the ./mapping directory should be executed first, and once this is finished, the scripts in the ./variant_calling directory can be run.

./testfiles contains example files (one chromosome from the Macaca mulatta reference genome, and a subset of reads extracted from a few publicly available primate genomes).

Documentation for the mapping and variant_calling processes are available in the respective README files located in each directory. 