#!/bin/bash
#PBS -N 3_processbam
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
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

# Set paths to input and output directories
export ALIGN=/user/2_Mapping/bowtie2_out
export RG=/user/3_Processing/read_groups
export DUP=/user/3_Processing/dedup
export SUM=/user/3_Processing/validate_bam
export REF=/user/2_Mapping/ref_data

# Setup array job
sed -n $PBS_ARRAYID"p" SAMPLES.txt > temp$PBS_ARRAYID
SAMPLE=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo Sample ID is $SAMPLE

# Add read groups and index file
gatk AddOrReplaceReadGroups \
-I $ALIGN/$SAMPLE'_'sorted.bam \
-O $RG/$SAMPLE'_'sorted_RG.bam \
--RGID=$PBS_ARRAYID \
--RGLB=lib$PBS_ARRAYID \
--RGPU=unit$PBS_ARRAYID \
--RGSM=$SAMPLE \
--RGPL=ILLUMINA \
--SORT_ORDER=coordinate \
--CREATE_INDEX=true \
--TMP_DIR=/user/3_Processing/TEMP

# Mark duplicates
gatk MarkDuplicates \
-I $RG/$SAMPLE'_'sorted_RG.bam \
-O $DUP/$SAMPLE'_'sorted_dedup_RG.bam \
-M $DUP/$SAMPLE'_'sorted_dedup_RG_metrics.txt \
--TMP_DIR=/user/3_Processing/TEMP

# Validate bam file
gatk ValidateSamFile \
-R $REF/Felis_catus_9_0.fa \
-I $DUP/$SAMPLE'_'sorted_dedup_RG.bam \
-O $SUM/$SAMPLE'_'validate.txt \
--TMP_DIR=/user/3_Processing/TEMP \
--MODE=SUMMARY

echo End time is `date`