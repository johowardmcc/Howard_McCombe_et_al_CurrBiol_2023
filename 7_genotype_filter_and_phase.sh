#!/bin/bash
#PBS -N 7_genotype
#PBS -l mem=40gb
#PBS -l nodes=1:ppn=10
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -t 1-20

echo "The Array ID is: $PBS_ARRAYID" 
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/gatk-4.0.8.1
module load apps/samtools-1.9
module load apps/tabix-0.2.6
module load apps/vcftools-0.1.17.2
module load apps/bcftools-1.8 

# Setup array job
sed -n $PBS_ARRAYID"p" CHROMOSOMES2.txt > temp$PBS_ARRAYID
CHROM=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo Chromosome $CHROM

# Set paths to input and output directories
export DB=/user/4_Genotyping/GenomicsDB
export REF=/user/2_Mapping/ref_data
export OUT=/user/4_Genotyping/GenotypeGVCF
export FILTER=/user/4_Genotyping/Filtered
export PHASE=/user/5_Phasing
export MAP=/user/5_Phasing/maps

# Genotype GVCFs
gatk GenotypeGVCFs \
	-R $REF/Felis_catus_9_0.fa \
	-V gendb://$DB/$PBS_ARRAYID \
	-O $OUT/$PBS_ARRAYID'.'vcf.gz \
	--TMP_DIR=/user/4_Genotyping/GenotypeGVCF/TEMP

# Process SNPs and indels separately (processing mixed-type variants with indels)
gatk SelectVariants \
    -V $OUT/$PBS_ARRAYID.vcf.gz \
    -select-type SNP \
    -O $OUT/$PBS_ARRAYID'_'snps.vcf.gz
	
gatk SelectVariants \
    -V $OUT/$PBS_ARRAYID.vcf.gz \
    -select-type INDEL \
    -select-type MIXED \
    -O $OUT/$PBS_ARRAYID'_'indels.vcf.gz
	
# GATK suggested hard filtering options
# SNPs
gatk VariantFiltration \
    -V $OUT/$PBS_ARRAYID'_'snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O $FILTER/$PBS_ARRAYID'_'snps_filtered.vcf.gz
	
# Indels
gatk VariantFiltration \
	-V $OUT/$PBS_ARRAYID'_'indels.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 200.0" --filter-name "FS200" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-O $FILTER/$PBS_ARRAYID'_'indels_filtered.vcf.gz
	
# Extract info about passed sites
gatk VariantsToTable \
    -V $FILTER/$PBS_ARRAYID'_'snps_filtered.vcf.gz \
	--show-filtered=false \
	-F QD -F QUAL -F SOR -F FS -F MQ -F MQRankSum -F ReadPosRankSum \
	-O $PLOTS/$PBS_ARRAYID'_'snps_PASS.txt

gatk VariantsToTable \
    -V $FILTER/$PBS_ARRAYID'_'indels_filtered.vcf.gz \
	--show-filtered=false \
	-F QD -F QUAL -F FS -F ReadPosRankSum \
	-O $PLOTS/$PBS_ARRAYID'_'indels_PASS.txt
		
# Extract info about all sites
gatk VariantsToTable \
    -V $FILTER/$PBS_ARRAYID'_'snps_filtered.vcf.gz \
	--show-filtered=true \
	-F QD -F QUAL -F SOR -F FS -F MQ -F MQRankSum -F ReadPosRankSum \
	-O $PLOTS/$PBS_ARRAYID'_'snps_ALL.txt

gatk VariantsToTable \
    -V $FILTER/$PBS_ARRAYID'_'indels_filtered.vcf.gz \
	--show-filtered=true \
	-F QD -F QUAL -F FS -F ReadPosRankSum \
	-O $PLOTS/$PBS_ARRAYID'_'indels_ALL.txt

# Remove GATK hard-filtered snps
bcftools view -i 'FILTER="PASS"' $FILTER/$PBS_ARRAYID'_'snps_filtered.vcf.gz -O z > $FILTER/$PBS_ARRAYID'_'GATK.vcf.gz
tabix -p vcf $FILTER/$PBS_ARRAYID'_'.GATK.vcf.gz

# Second round of SNP filtering
vcftools --gzvcf $FILTER/$PBS_ARRAYID'_'.GATK.vcf.gz \
	--keep FINAL_SAMPLES.txt \
	--max-missing 1 \
	--minQ 49 \
	--maxDP 2000 \
	--mac 3 \
	--min-alleles 2 --max-alleles 2 \
	--recode --recode-INFO-all --stdout | bgzip -c > $FILTER/$PBS_ARRAYID'_'.FILTER2.vcf.gz

# Phase
java -jar beagle.jar gt=$FILTER/$PBS_ARRAYID'_'.FILTER2.vcf.gz map=$MAP/${CHROM}.MAP.txt burnin=50 iterations=100 out=$PHASE/${CHROM}.PHASED nthreads=8

echo End time is `date`
