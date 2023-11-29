require(MOSAIC)

args = commandArgs(trailingOnly=TRUE)

# Run MOSAIC
mosaic.result=run_mosaic("scot","/user/analysis/mosaic/INPUT/",1:19,A=2,NUMI=36)

# Generate bootstrap confidence intervals
# Load results
load("MOSAIC_RESULTS/scot_2way_1-36_1-19_130_60_0.99_58.RData")
load("MOSAIC_RESULTS/localanc_scot_2way_1-36_1-19_130_60_0.99_58.RData")

# Bootstrap (default 100 iterations)
bfit=bootstrap_chromosomes_coanc_curves(acoancs, dr, localanc, alpha)

x1<-numeric()
x2<-numeric()
x3<-numeric()
n=length(bfit$boot.gens)

for (j in 1:n){
	x<-data.frame(bfit$boot.gens[j])
	x1[j]<-x[1,1]
	x2[j]<-x[1,2]
	x3[j]<-x[2,2]
	}

bootstrap<-data.frame("1_1"=x1, "1_2"=x2, "2_2"=x2)
write.table(bootstrap, "bootstrap2.txt",quote=F, col.names=T, row.names=F)
