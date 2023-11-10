require(MOSAIC)
library(vcfR)

load("masking_data/scot_2way_1-36_1-19_130_60_0.99_58.RData")
load("masking_data/localanc_scot_2way_1-36_1-19_130_60_0.99_58.RData")

unphased_localanc=phase_localanc(localanc,final.flips) 
local_pos=grid_to_pos(unphased_localanc,"masking_data/INPUT2/",g.loc,chrnos)


for (chrom in 1:18){
	wildcat_data<-read.vcfR(paste("chr",chrom,".recode.vcf", sep=""))
	array1=local_pos[[chrom]]
	nsnps<-dim(array1)[3]
	nindv<-dim(array1)[2]/2 

	for (ind in 1:nindv){
		mask1<-numeric()	
		mask2<-numeric()
		hap1=ind*2-1
		hap2=hap1+1	
	
		# This checks the probability of wildcat ancestry 
		# at both haplotypes of the individual, at each SNP. If either haplotype 
		# has less than a 80% chance of coming from a wildcat we want to remove this SNP,
		# so we set the 'masking vector' to 0.
	
		for (snp in 1:nsnps) {
			if (array1[2, hap1, snp] > 0.8) {
				mask1[snp]<- 1 
				} else {
				mask1[snp]<- 0
				}
		
			if (array1[2, hap2, snp] > 0.8) {
				mask2[snp]<- 1 
				} else {
				mask2[snp]<- 0
				}
			}
	
		# Replace 'domestic' sites with missing data (".") in the vcf using the information in 'mask'
	
		gvec1 = substr(wildcat_data@gt[,ind+8],1,1) # all the sites from 1st copy of the chromosome
		gvec1 = ifelse(mask1,gvec1,".")
		gvec2 = substr(wildcat_data@gt[,ind+8],3,3) # all the sites from the 2nd copy
		gvec2 = ifelse(mask2,gvec2,".")
		
		# Join everything together to create a genotype per site
		gvec_final = paste(gvec1,"|",gvec2,sep="") # combine the two copies with the separator "|"

		# Replace original genotypes with masked set (gvec_final)
		substr(wildcat_data@gt[,ind+8],1,3) <- gvec_final

	}
	
	write.vcf(wildcat_data, file=paste("chr",chrom,".masked.vcf.gz",sep=""))
}