#!/bin/bash
#PBS -N pca
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=2
#PBS -l walltime=10:00:00
#PBS -j oe

echo "The Array ID is: $PBS_ARRAYID"
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs 
module load apps/eigensoft-7.2.1
module load apps/samtools-1.9
module load apps/tabix-0.2.6
module load apps/vcftools-0.1.17.2
module load apps/bcftools-1.8
module load apps/plink-1.90

# Set paths to input and output directories
export WD=/user/analysis
export MODERN=/user/analysis/wildcats_final_phased.vcf.gz
export HISTORIC=/user/analysis/historic_wgs_final.vcf.gz
export ANCIENT=/user/analysis/ancient_silvestris.vcf.gz

cd $WD

bcftools merge -m snps -O z --threads 2 -o temp.vcf.gz $MODERN $HISTORIC $ANCIENT 
tabix temp.vcf.gz

vcftools --gzvcf temp.vcf.gz \
	--positions common_sites_final.txt \
	--min-alleles 2 --max-alleles 2 \
	--recode --stdout | bgzip -c > eigensoft.vcf.gz
tabix eigensoft.vcf.gz

bcftools query -l eigensoft.vcf.gz > eigensoft_samples.txt

rm temp.vcf.gz*

plink --vcf eigensoft.vcf.gz --chr-set 18 --keep-allele-order --set-missing-var-ids @:# --bp-space 1000 --pheno pheno.file --recode --make-bed --out eigensoft

# PCA with eigensoft
convertf -p parfile.convert

smartpca -p parfile.smartpca
