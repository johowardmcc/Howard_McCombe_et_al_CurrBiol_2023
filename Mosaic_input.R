# Aug2021

args = commandArgs(trailingOnly=TRUE)
chrom=(args[1])
path1="/newhome/jh17379/CHAPTER5/6_MOSAIC/INPUT2/"

if (chrom==1){

        # Population name file
        dom<-read.table(paste(path1,chrom,".DOM.impute.hap.indv",sep=""))
        scot<-read.table(paste(path1,chrom,".SCOT.impute.hap.indv",sep=""))
        wildcat<-read.table(paste(path1,chrom,".WILDCAT.impute.hap.indv",sep=""))
	otherdom<-read.table(paste(path1,chrom,".OTHERDOM.impute.hap.indv",sep=""))

        pop.file<-data.frame(pop=c(rep("scot",dim(scot)[1]),
		rep("dom",dim(dom)[1]),
		rep("otherdom",dim(otherdom)[1]),
		rep("wildcat",dim(wildcat)[1])),
                indv=rbind(scot,dom,otherdom,wildcat))
        write.table(pop.file,paste(path1,"sample.names",sep=""),row.names=F,col.names=F,quote=F)
        }

# Make SNP file
data<-read.table(paste(path1,chrom,".SCOT.impute.legend",sep=""),header=T)

snp.file<-data.frame(data$ID,
        rep(chrom,dim(data)[1]),
        rep("?",dim(data)[1]),
        data$pos,
        data$allele0,
        data$allele1)

write.table(snp.file,paste(path1,"snpfile.",chrom,sep=""),quote=F,row.names=F,col.names=F)
