#!/bin/bash
#PBS -N 2_alignment
#PBS -l mem=32gb
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -t 1-45

echo "The Array ID is: $PBS_ARRAYID"
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/bowtie-2.2.9
module load apps/samtools-1.9
module load apps/fastqc-0.11.2

# set paths to input, output and reference directories
export REF=/user/2_Mapping/ref_data
export TRIM=/user/1_QC/trim_out
export OUT=/user/2_Mapping/bowtie2_out

# Setup array job
sed -n $PBS_ARRAYID"p" SAMPLES.txt > temp$PBS_ARRAYID
SAMPLE=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo Sample ID is $SAMPLE

# align trimmed reads with bowtie2
bowtie2 -1 $TRIM/$SAMPLE'_'trimmed_1P.fq.gz -2 $TRIM/$SAMPLE'_'trimmed_2P.fq.gz -U $TRIM/$SAMPLE'_'trimmed_1U.fq.gz,$TRIM/$SAMPLE'_'trimmed_2U.fq.gz -x $REF/Felis_catus -p 8 -t -q | samtools sort -m 32G -O BAM -o $OUT/$SAMPLE'_'sorted.bam

# Set path to QC directory
export QC1=/user/2_Mapping/samtools_flagstats
export QC2=/user/2_Mapping/fastqc

# Get alignment stats
cd $QC1
samtools flagstat $OUT/$SAMPLE'_'sorted.bam 

cd ../..

# Fastqc sorted bam 
fastqc $OUT/$SAMPLE'_'sorted.bam -o $QC2

echo End time is `date`

