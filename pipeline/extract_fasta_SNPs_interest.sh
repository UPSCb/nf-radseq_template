#!/bin/bash
#SBATCH --account=naiss2024-22-1261
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --time=1:00:00

grep "#" analysis/stacks/populations.snps.vcf > analysis/interest_SNPs.vcf
grep -v "INFO" analysis/interest_SNPs.tsv | cut -f 27,28,29,30,31,32,33,34,35,36,37 --complement >> analysis/interest_SNPs.vcf


module load samtools
grep -v "#" analysis/stacks/populations.loci.fa > analysis/stacks/populations.loci.fasta
rm analysis/interest_SNPs.fasta

for i in $(cut -f 37 analysis/interest_SNPs.tsv | grep -v coordinates)
do samtools faidx analysis/stacks/populations.loci.fasta "$i" >> analysis/interest_SNPs.fasta
done


/cfs/klemming/home/e/edpi01/miniforge3/envs/biopandas/bin/python src/python/Assemble_final_output.py analysis/interest_SNPs.tsv \
analysis/interest_SNPs.fasta analysis/SNPs_interest_WithSequences.csv 5
