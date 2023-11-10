#!/bin/bash
#PBS -N 8_concatdata
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00
#PBS -j oe

echo "The Array ID is: $PBS_ARRAYID" 
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/samtools-1.9
module load apps/tabix-0.2.6
module load apps/vcftools-0.1.17.2
module load apps/bcftools-1.8 

# Set paths to input and output directories
export PHASE=/user/5_Phasing
export FINAL=/user/analyses

# concat data from all chromosomes
bcftools concat --file-list phased.vcf.txt -O z > temp.vcf.gz
tabix -p vcf temp.vcf.gz

# change chromosome names to numbers
bcftools annotate --rename-chr name2number.txt -O z -o wildcats_final_phased.vcf.gz temp.vcf.gz
tabix -p vcf wildcats_final_phased.vcf.gz

# make plink and thin
plink --vcf wildcats_final_phased.vcf.gz \
	--chr-set 18 \
	--keep-allele-order \
	--set-missing-var-ids @:# \
	--bp-space 2000 \
	--make-bed --out wildcats_final_phased_thinned

rm temp.vcf.*