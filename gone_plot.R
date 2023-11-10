catjack = function(jackvals,pval=0.95)
#jackvals is the vector of estimates from a leave-one-out procedure
#
#output: vector with (jackknife mean,jackknife se, jackkknife confidence interval)
#i.e. if the jackknife mean is mu and confidence interval is ci 
#then the summary is mu +/- ci

{
	n = length(jackvals)
	jackmean = mean(jackvals)
	sejack <- sqrt(((n - 1)/n) * sum((jackvals - mean(jackvals))^2))
	crit.val = qnorm(1 - (1-pval)*0.5) # note: 2 tailed
	cijack = crit.val*sejack
	return(c(jackmean,sejack,cijack))
	
} 

# Jackknife bias calculated using: (n - 1) * (mean(leave_one_out_estims) - estim_all), including information from a GONE run using all chromosomes


## WILDCAT ##
wild.data<-read.table("WILDCAT/Output_Ne_1",header=T)
wild.data<-wild.data[1:150,]
colnames(wild.data)<-c("GEN","RUN1")

for (j in 2:18){
	run<-read.table(paste("WILDCAT/Output_Ne_",j,sep=""), header=T)
	run<-run[1:150,]
	wild.data[,j+1]<-run[,2]
	}
	
mean.ne<-numeric()
upr<-numeric()
lwr<-numeric()
for (j in 1:150){
	row<-as.numeric(wild.data[j,])
	row<-row[-1]
	stats<-catjack(row)
	mean.ne[j]<-stats[1]
	upr[j]<-stats[1]+stats[3]
	lwr[j]<-stats[1]-stats[3]
	}

wild.data$WILD.MEAN<-mean.ne
wild.data$WILD.LWR<-lwr
wild.data$WILD.UPR<-upr
colnames(wild.data)<-c("GEN","WILD.RUN1","WILD.RUN2","WILD.RUN3","WILD.RUN4","WILD.RUN5","WILD.RUN6","WILD.RUN7","WILD.RUN8","WILD.RUN9","WILD.RUN10","WILD.RUN11","WILD.RUN12","WILD.RUN13","WILD.RUN14","WILD.RUN15","WILD.RUN16","WILD.RUN17","WILD.RUN18","WILD.MEAN","WILD.LWR","WILD.UPR")

## DOMESTIC ##
dom.data<-read.table("DOMESTIC/Output_Ne_1",header=T)
dom.data<-dom.data[1:150,]
colnames(dom.data)<-c("GEN","RUN1")

for (j in 2:18){
	run<-read.table(paste("DOMESTIC/Output_Ne_",j,sep=""), header=T)
	run<-run[1:150,]
	dom.data[,j+1]<-run[,2]
	}
	
mean.ne<-numeric()
upr<-numeric()
lwr<-numeric()
for (j in 1:150){
	row<-as.numeric(dom.data[j,])
	row<-row[-1]
	stats<-catjack(row)
	mean.ne[j]<-stats[1]
	upr[j]<-stats[1]+stats[3]
	lwr[j]<-stats[1]-stats[3]
	}

dom.data$DOM.MEAN<-mean.ne
dom.data$DOM.LWR<-lwr
dom.data$DOM.UPR<-upr

colnames(dom.data)<-c("GEN","DOM.RUN1","DOM.RUN2","DOM.RUN3","DOM.RUN4","DOM.RUN5","DOM.RUN6","DOM.RUN7","DOM.RUN8","DOM.RUN9","DOM.RUN10","DOM.RUN11","DOM.RUN12","DOM.RUN13","DOM.RUN14","DOM.RUN15","DOM.RUN16","DOM.RUN17","DOM.RUN18","DOM.MEAN","DOM.LWR","DOM.UPR")
dom.data<-dom.data[,-1]

## PLOT ##
data<-cbind(wild.data,dom.data)

for (j in 1:150){
	if (data$WILD.LWR[j] < 0) {
		data$WILD.LWR[j]<-0
		}
	if (data$DOM.LWR[j] < 0) {
		data$DOM.LWR[j] <- 0
		}
	}

x<-ggplot(data)+
  geom_line(aes(GEN, WILD.MEAN, col="ancestry1"), size=1.2)+
	geom_line(aes(GEN, DOM.MEAN, col="ancestry2"),size=1.2)+
	geom_ribbon(aes(x=GEN, ymin=WILD.LWR, ymax=WILD.UPR),fill="#56B4E9",alpha=0.2)+
	geom_ribbon(aes(x=GEN, ymin=DOM.LWR, ymax=DOM.UPR),fill="#E69F00",alpha=0.2)+
	scale_color_manual(name="ancestry", labels=c("Wildcat ancestry","Domestic ancestry"), values=c("#56B4E9", "#E69F00"))+
	xlab("Generations before present")+
	ylab("Ne")+
	theme_classic()+
	theme(legend.title=element_blank(),legend.direction="horizontal",legend.position = c(0.5, 0.97))+
  guides(color=guide_legend(override.aes=list(fill=NA)))

## ADMIXED DATA ##
admix.data<-read.table("ADMIXED/Output_Ne_1",header=T)
admix.data<-admix.data[1:150,]
colnames(admix.data)<-c("GEN","RUN1")

for (j in 2:18){
	run<-read.table(paste("ADMIXED/Output_Ne_",j,sep=""), header=T)
	run<-run[1:150,]
	admix.data[,j+1]<-run[,2]
	}
	
mean.ne<-numeric()
upr<-numeric()
lwr<-numeric()
for (j in 1:150){
	row<-as.numeric(admix.data[j,])
	row<-row[-1]
	stats<-catjack(row)
	mean.ne[j]<-stats[1]
	upr[j]<-stats[1]+stats[3]
	lwr[j]<-stats[1]-stats[3]
	}

admix.data$ADMIX.MEAN<-mean.ne
admix.data$ADMIX.LWR<-lwr
admix.data$ADMIX.UPR<-upr

## PLOT ##
for (j in 1:150){
	if (admix.data$ADMIX.LWR[j] < 0) {
		admix.data$ADMIX.LWR[j]<-0
		}
	}

z<-ggplot(admix.data)+	
	geom_line(aes(GEN, ADMIX.MEAN, col="Admixed"), size=1.2)+
	geom_ribbon(aes(x=GEN, ymin=ADMIX.LWR, ymax=ADMIX.UPR),fill="#289C61",alpha=0.2)+
	scale_color_manual(name="Admixed", labels="Complete (admixed) data", values="#289C61")+
	xlab("Generations before present")+
	ylab("Ne")+
	theme_classic()+
	theme(legend.title=element_blank(),legend.position = c(0.3, 0.97))

grid.arrange(x,z,ncol=2)
dev.off()