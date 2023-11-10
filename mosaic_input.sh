#PBS -N mosaic1
#PBS -l mem=5gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -t 1-19

echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID

# Load programs
module load apps/vcftools-0.1.17.0
module load apps/samtools-1.9.1
module load apps/tabix-0.2.6
module load languages/R-3.5.1-ATLAS-gcc-6.1

# Set paths to relevant directories
export WD=/user/analysis
export VCF=/user/analysis/wildcats_final_phased_thinned.vcf.gz

cd $WD

# Make genofile

vcftools --gzvcf $VCF --keep scotland.txt --chr $PBS_ARRAYID --recode --stdout | bgzip -c > ${PBS_ARRAYID}.SCOT.vcf.gz
vcftools --gzvcf $VCF --keep domestic.txt --chr $PBS_ARRAYID --recode --stdout | bgzip -c > ${PBS_ARRAYID}.DOM.vcf.gz
vcftools --gzvcf $VCF --keep otherdom.txt --chr $PBS_ARRAYID --recode --stdout | bgzip -c > ${PBS_ARRAYID}.OTHERDOM.vcf.gz
vcftools --gzvcf $VCF --keep wildcat.txt --chr $PBS_ARRAYID --recode --stdout | bgzip -c > ${PBS_ARRAYID}.WILDCAT.vcf.gz

vcftools --gzvcf ${PBS_ARRAYID}.SCOT.vcf.gz --IMPUTE --out ${PBS_ARRAYID}.SCOT
vcftools --gzvcf ${PBS_ARRAYID}.DOM.vcf.gz --IMPUTE --out ${PBS_ARRAYID}.DOM
vcftools --gzvcf ${PBS_ARRAYID}.OTHERDOM.vcf.gz --IMPUTE --out ${PBS_ARRAYID}.OTHERDOM
vcftools --gzvcf ${PBS_ARRAYID}.WILDCAT.vcf.gz --IMPUTE --max-missing 0 --out ${PBS_ARRAYID}.WILDCAT

cat ${PBS_ARRAYID}.SCOT.impute.hap | tr -d "[:blank:]" > scotgenofile.${PBS_ARRAYID}
cat ${PBS_ARRAYID}.DOM.impute.hap | tr -d "[:blank:]" > domgenofile.${PBS_ARRAYID}
cat ${PBS_ARRAYID}.OTHERDOM.impute.hap | tr -d "[:blank:]" > otherdomgenofile.${PBS_ARRAYID}
cat ${PBS_ARRAYID}.WILDCAT.impute.hap | tr -d "[:blank:]" > wildcatgenofile.${PBS_ARRAYID}

Rscript Mosaic_input.R $PBS_ARRAYID

rm ${PBS_ARRAYID}.SCOT.*
rm ${PBS_ARRAYID}.DOM.*
rm ${PBS_ARRAYID}.OTHERDOM.*
rm ${PBS_ARRAYID}.WILDCAT.*
