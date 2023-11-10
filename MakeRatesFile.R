args = commandArgs(trailingOnly=TRUE)
chrom=(args[1])

path="/newhome/jh17379/CHAPTER5/6_MOSAIC/INPUT2/"

map<-read.table(paste(path,chrom,".MOSAIC.txt",sep=""),header=T)

map$cM=map[,2]*100

end=dim(map)[1]
map[end,3]<-map[end-1,3]

map<-data.frame(map[,1],map[,3])

write.table(t(map),paste(path,"rates.",chrom,sep=""),quote=F,row.names=F,col.names=F)
