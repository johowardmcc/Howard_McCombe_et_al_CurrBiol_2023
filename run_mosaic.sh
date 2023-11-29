#!/bin/bash
#PBS -N Mosaic
#PBS -l mem=30gb
#PBS -l nodes=1:ppn=8
#PBS -l walltime=50:00:00
#PBS -j oe

# OPTIONAL: record some potentially useful details about the job:
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID

#load your program if it is installed globally
module load languages/R-3.5.1-ATLAS-gcc-6.1

# set a local temporary folder
export TMPDIR=$HOME/.local

# Set paths to relevant directories
export WD=/user/analysis/mosaic

cd $WD

# Run mosaic
Rscript RunMosaic.R

echo End time is `date`
