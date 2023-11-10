#!/bin/bash
#PBS -N mosaic2
#PBS -l mem=5gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -t 1-19

echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID

# Load programs
module load apps/vcftools-0.1.17.0
module load apps/plink-1.90
module load languages/R-3.5.1-ATLAS-gcc-6.1

# Set paths to relevant directories
export WD=/user/analysis
export VCF=/user/analysis/wildcats_final_phased_thinned.vcf.gz
export HAPMAP=/user/analysis/HAPMAP

cd $WD

# Use FS to make map
vcftools --gzvcf $VCF --chr $PBS_ARRAYID --plink --temp TEMP --out $PBS_ARRAYID
plink --file $PBS_ARRAYID --keep-allele-order --chr-set 18 --recode 12 --out $PBS_ARRAYID

# executables from finestructure (
plink2chromopainter.pl -p=${PBS_ARRAYID}.ped -m=${PBS_ARRAYID}.map -o=${PBS_ARRAYID}.phase

convertrecfile.pl -U c -t 1,3 -D CDF -c ${PBS_ARRAYID}.phase $HAPMAP/${PBS_ARRAYID}.HapMap.txt ${PBS_ARRAYID}.MOSAIC.txt

# Convert map file
Rscript MakeRatesFile.R $PBS_ARRAYID

## N.B. Output rates files need to be edited to include the number of sites included as the first line! ##

rm ${PBS_ARRAYID}.log
rm ${PBS_ARRAYID}.map
rm ${PBS_ARRAYID}.nosex
rm ${PBS_ARRAYID}.ped
rm ${PBS_ARRAYID}.phase

echo End time is `date`
