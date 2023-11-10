require(MOSAIC)
library(vcfR)

load("masking_data/scot_2way_1-36_1-19_130_60_0.99_58.RData")
load("masking_data/localanc_scot_2way_1-36_1-19_130_60_0.99_58.RData")

unphased_localanc=phase_localanc(localanc,final.flips) 
local_pos=grid_to_pos(unphased_localanc,"/user/analysis/",g.loc,chrnos)


for (chrom in 1:18){
	wildcat_data<-read.vcfR(paste("masked_chr/chr",chrom,".recode.vcf", sep=""))
	array1=local_pos[[chrom]]
	nsnps<-dim(array1)[3]
	nindv<-dim(array1)[2]/2 
	
	ID<-colnames(wildcat_data@gt)
	col.names<-vector()
	col.names[1]<-"FORMAT"
	
	wildcat_hap<-data.frame(FORMAT=c(rep("GT",nsnps)))
	eur_hap<-data.frame(FORMAT=c(rep("GT",nsnps)))
	eur.names<-vector()
	dom_hap<-data.frame(FORMAT=c(rep("GT",nsnps)))
	dom.names<-vector()
	
	for (ind in 1:nindv){
		ind.id<-ID[ind+8]
		print(ind.id)
		mask1<-numeric()	
		mask2<-numeric()
		hap1=ind*2-1
		hap2=hap1+1	
	
		# This loop within a loop checks the probability of wildcat ancestry 
		# at both haplotypes of the individual, at each SNP. If either haplotype 
		# has less than a 80% chance of coming from a wildcat we want to remove this SNP,
		# so we set the 'masking vector' to 0.
		# At the end of this loop we should have a vector called 'mask' which is the same 
		# length as the number of SNPs on the chromosome, with 0s for the sites we want to
		# remove and 1s for the ones we want to keep
	
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
		
		col.names[ind*2]<-paste(ind.id,"_A",sep="")
		col.names[(ind*2)+1]<-paste(ind.id,"_B",sep="")
		
		#gvec1_final<-paste(gvec1,"|.",sep="")
		#gvec2_final<-paste(".|",gvec2,sep="")
		
		wildcat_hap[,ind*2]<-gvec1
		wildcat_hap[(ind*2)+1]<-gvec2
		
	}
	
	for (ind in 1:7){
	
		gvec1<-substr(wildcat_data@gt[,ind+1],1,1)
		gvec2<-substr(wildcat_data@gt[,ind+1],3,3)
		
		eur_hap[,(ind*2)]<-gvec1
		eur_hap[,(ind*2)+1]<-gvec2
		
		ind.id<-ID[ind+1]
		print(ind.id)
		eur.names[(ind*2)-1]<-paste(ind.id,"_A",sep="")
		eur.names[(ind*2)]<-paste(ind.id,"_B",sep="")
		
	}
	
	for (ind in 1:22){
	
		gvec1<-substr(wildcat_data@gt[,ind+44],1,1)
		gvec2<-substr(wildcat_data@gt[,ind+44],3,3)
		
		dom_hap[,(ind*2)]<-gvec1
		dom_hap[,(ind*2)+1]<-gvec2
		
		ind.id<-ID[ind+44]
		print(ind.id)
		dom.names[(ind*2)-1]<-paste(ind.id,"_A",sep="")
		dom.names[(ind*2)]<-paste(ind.id,"_B",sep="")
		
	}

	eur_hap<-eur_hap[,-1]
	colnames(eur_hap)<-eur.names
	dom_hap<-dom_hap[,-1]
	colnames(dom_hap)<-dom.names
	colnames(wildcat_hap)<-col.names
	all_hap<-cbind(wildcat_hap,eur_hap,dom_hap)
	all_hap<-as.matrix(all_hap)
	wildcat_data@gt<-all_hap
	
	write.vcf(wildcat_data, file=paste("haploid/chr",chrom,".masked.vcf.gz",sep=""))
	
	print(paste("Finished chromosome ", chrom, sep="")) 
}


