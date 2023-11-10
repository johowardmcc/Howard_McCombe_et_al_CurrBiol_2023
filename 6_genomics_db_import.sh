#!/bin/bash
#PBS -N 6_genomicsdb
#PBS -l mem=40gb
#PBS -l nodes=1:ppn=10
#PBS -l walltime=150:00:00
#PBS -j oe
#PBS -t 1-20 # number of chromosomes + mtDNA

echo "The Array ID is: $PBS_ARRAYID" 
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load programs
module load apps/gatk-4.0.8.1

# Setup array job
sed -n $PBS_ARRAYID"p" CHROMOSOMES.txt > temp$PBS_ARRAYID
CHROM=$(awk '{print $1}' temp$PBS_ARRAYID)
rm temp$PBS_ARRAYID

echo CHROM ID is $CHROM

# Set paths to input and output directories
export OUT=/user/4_Genotyping/GenomicsDB

gatk GenomicsDBImport \
	--sample-name-map GVCF_PATH.txt \
	--genomicsdb-workspace-path $OUT/$PBS_ARRAYID \
	--TMP_DIR /user/4_Genotyping/TEMP \
	-L $CHROM \
	--batch-size 20 \
	--reader-threads 10
	
echo End time is `date`

