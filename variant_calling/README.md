# Haplotype Caller and Joint Genotyping Pipeline    

In a first step, HaplotypeCaller is run per sample. A second step, which is executed automatically, combines all sample-specific haplotype calls using CombineGVCFs, and subsequently genotypes them jointly using GenotypeGVCFs.

The output of the pipeline will be jointly genotyped vcf files, one per interval/chunk as specified in the start_genotyping.sh script. These can then be filtered depending on use case and combined into chromosome or genome level vcf files with e.g. bcftools concat.

### Dependencies

The script uses the following software:
- samtools/1.14
- GATK/4.3.0.0

These are available as preinstalled modules on the cluster for which this pipeline was developed (Uppmax).

### Running the pipeline

1) Set up a file containing the scaffolds/chromosomes upon which to call variants. The pipeline will be run separately for each scaffold/chromosome, and to increase speed/parallelization I've typically split up large scaffolds/chromosomes into smaller chunks (~20 Mb). Note that if the reference contains many tiny, unplaced contigs/scaffolds, it might be a good idea to exclude them already at this step, since many jobs will be started per sample otherwise. The fils should have at least two tab-separated columns: scaffold and length. A faidx file works too, everything beyond the first two cols will be ignored.

2) Modify the parameters in the start_genotyping.sh script, such that the paths to the reference genome, directory with mapped bam files etc. are correct.

3) Start the pipeline simply by running the start_genotyping.sh script:
```
bash start_genotyping.sh
```
This will spawn array HaplotypeCaller jobs for each interval, with one task per sample. Once all HaplotypeCaller jobs are finished, the second step with CombineGVCFs and GenotypeGVCFs will be started automatically.

