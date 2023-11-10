#!/bin/bash
#PBS -N 5_haplotypecaller
#PBS -l mem=40gb
#PBS -l nodes=1:ppn=10
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -t 1-45

echo "The Array ID is: $PBS_ARRAYID" 
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/gatk-4.0.8.1
module load apps/samtools-1.9.1

# Setup array job
sed -n $PBS_ARRAYID"p" SAMPLES.txt > temp$PBS_ARRAYID
SAMPLE=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo Sample ID is $SAMPLE

# Set paths to input and output directories
export FINAL=/user/3_Processing/fix_bam
export REF=/user/2_Mapping/ref_data
export OUT=/user/4_Genotyping/GVCF

# Make sure all BAMs files are indexed
samtools index $FINAL/$SAMPLE'_'sorted_dedup_RG_fixed_mate.bam 

gatk HaplotypeCaller \
	-R $REF/Felis_catus_9_0.fa \
	-I $FINAL/$SAMPLE'_'sorted_dedup_RG_fixed_mate.bam \
	-ERC GVCF \
	-O $OUT/$SAMPLE'_'final.g.vcf \
	--TMP_DIR=/user/4_Genotyping/TEMP

echo End time is `date`