#!/bin/bash
#PBS -N 1_QC
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=8
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -t 1-45

echo "The Array ID is: $PBS_ARRAYID"
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/fastqc-0.11.2
module load languages/java-jdk-1.7.0-40

# Setup array job
sed -n $PBS_ARRAYID"p" SAMPLES.txt > temp$PBS_ARRAYID
SAMPLE=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo Sample ID is $SAMPLE

# Concatenate reads
cat $SAMPLE/*_1.fq.gz > 1_QC/raw_data/$SAMPLE_1.fq.gz
cat $SAMPLE/*_2.fq.gz > 1_QC/raw_data/$SAMPLE_2.fq.gz

# FastQC where raw data is stored
cd 1_QC/raw_data

fastqc $SAMPLE_1.fq.gz $SAMPLE_2.fq.gz -o ../fastqc_out

# Trimmomatic where raw data is stored
# trimmomatic-0.39.jar must also be in this folder
export OUT=1_QC/trim_out

java -Xmx32g -jar ./trimmomatic-0.39.jar PE -threads 4 -phred33 \
	$SAMPLE_1.fq.gz $SAMPLE_2.fq.gz \
	$OUT/$SAMPLE_trimmed_1P.fq.gz  $OUT/$SAMPLE_trimmed_1U.fq.gz  $OUT/$SAMPLE_trimmed_2P.fq.gz  $OUT/$SAMPLE_trimmed_2U.fq.gz \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo End time is `date`

