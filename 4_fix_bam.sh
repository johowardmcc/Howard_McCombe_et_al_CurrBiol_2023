#!/bin/bash
#PBS -N 4_fixbam
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -t 1-45 # only use per sample as required

echo "The Array ID is: $PBS_ARRAYID" # this one is for when you are running an array job
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/gatk-4.0.8.1

# Setup array job
sed -n $PBS_ARRAYID"p" SAMPLES.txt > temp$PBS_ARRAYID
SAMPLE=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo Sample ID is $SAMPLE

# Set paths to input and output directories
export REF=/user/2_Mapping/ref_data
export DUP=/user/3_Processing/dedup
export OUT=/user/3_Processing/fix_bam
export SUM=/user/3_Processing/validate_bam

# Fix NM tag
gatk SetNmMdAndUqTags \
	-R=$REF/Felis_catus_9_0.fa \
	-I=$DUP/$SAMPLE'_'sorted_dedup_RG.bam \
	-O=$OUT/$SAMPLE'_'sorted_dedup_RG_fixed_NM.bam \
	--TMP_DIR=$OUT/TEMP

# Fix mate information
gatk FixMateInformation \
	-I=$OUT/$SAMPLE'_'sorted_dedup_RG_fixed_NM.bam \
	-O=$OUT/$SAMPLE'_'sorted_dedup_RG_fixed_mate.bam \
	--ADD_MATE_CIGAR=true \
	--SORT_ORDER=coordinate \
	--TMP_DIR=$OUT/TEMP

# Keep final fixed BAM	
rm $OUT/$SAMPLE'_'sorted_dedup_RG_fixed_NM.bam 

# Validate bam file again
gatk ValidateSamFile \
	-R $REF/Felis_catus_9_0.fa \
	-I $OUT/$SAMPLE'_'sorted_dedup_RG_fixed_mate.bam \
	-O $SUM/$SAMPLE'_'validate2.txt \
	--CREATE_INDEX=true \
	--TMP_DIR=$SUM/TEMP \
	--MODE=SUMMARY


echo End time is `date`