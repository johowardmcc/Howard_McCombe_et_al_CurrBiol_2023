#!/bin/bash
#PBS -N gone
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=2
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -t 1-18

echo "The Array ID is: $PBS_ARRAYID"
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This job runs on the following machines:

# Load your program
module load apps/plink/1.90 

# Record sample ID
echo Chromosome number $SLURM_ARRAY_TASK_ID

## GONE is run with a leave-one-out jackknife approach.  It take as input the masked, haploid data (using information from mosaic) for
## ancestry-specific analyses, or the full (admixed) data for Scottish wildcats.  It uses the default GONE settings in the INPUT_PARAMETERS_FILE
## except with KNOWN PHASE (1) for the admixed data and PSEUDOHAPLOID (0) for the ancestry-specific analyses.
## Reombination information was added to the .map files using same rate files as those generated for MOSAIC.

# Domestic
cd DOMESTIC
mkdir $SLURM_ARRAY_TASK_ID
cd $SLURM_ARRAY_TASK_ID
cp ../INPUT_PARAMETERS_FILE .
cp ../script_GONE.sh .
cp -r ../PROGRAMMES/ .

plink --file ../domestic.haploid --chr-set 18 --keep-allele-order --not-chr $SLURM_ARRAY_TASK_ID --recode --out $SLURM_ARRAY_TASK_ID

bash script_GONE.sh $SLURM_ARRAY_TASK_ID

# Wildcat
cd ../../WILDCAT
mkdir $SLURM_ARRAY_TASK_ID
cd $SLURM_ARRAY_TASK_ID
cp ../INPUT_PARAMETERS_FILE .
cp ../script_GONE.sh .
cp -r ../PROGRAMMES/ .

plink --file ../wildcat.haploid --chr-set 18 --keep-allele-order --not-chr $SLURM_ARRAY_TASK_ID --recode --out $SLURM_ARRAY_TASK_ID

bash script_GONE.sh $SLURM_ARRAY_TASK_ID

# Admixed
cd ../../ADMIXED
mkdir $SLURM_ARRAY_TASK_ID
cd $SLURM_ARRAY_TASK_ID
cp ../INPUT_PARAMETERS_FILE .
cp ../script_GONE.sh .
cp -r ../PROGRAMMES/ .

plink --file ../admix.scotland --chr-set 18 --keep-allele-order --not-chr $SLURM_ARRAY_TASK_ID --recode --out $SLURM_ARRAY_TASK_ID

bash script_GONE.sh $SLURM_ARRAY_TASK_ID


