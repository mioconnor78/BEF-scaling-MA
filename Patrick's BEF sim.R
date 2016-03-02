nCom<-100

initial<-round(rnorm(nCom,mean = 15,sd=2.5))

#dist<-(rbeta(n = nCom,shape1 = 5,shape2 = 2)*30)
#change<-round(dist-mean(dist))

change<-round(rnorm(nCom,mean = 0,sd = 2))

final<-initial+change
final[final<0]<-0

a<-4.33
b<-0.26

bmass1<-exp(a+b*log(initial))
bmass2<-exp(a+b*log(final))

bmass_change<-(bmass2-bmass1)/bmass1



pdf("Biomass change sim.pdf",width = 10,height = 6)
par(mfrow=c(1,2))
hist(change,xlab="Richness change",col="dark green", main="",breaks=seq(min(change)-1,max(change),by=1))
hist(bmass_change,col="grey35",right = T,breaks=seq(-round(max(abs(bmass_change)),2)-0.02,round(max(abs(bmass_change)),2)+0.02,by=0.02),main="",xlab="Biomass change")
dev.off()

sim_change<-function(nCom=1000,N_mean=10,N_sd=3,change=-0.2,a=4.33,b=0.26){
  initial<-round(rnorm(nCom,mean = N_mean,sd = N_sd))
  initial[initial<0]<-0
  final<-round(rnorm(nCom,mean = N_mean*(1+change),sd = N_sd))
  final[final<0]<-0
  bmass1<-exp(a+b*log(initial))
  bmass2<-exp(a+b*log(final))
  bmass_change<-(bmass2-bmass1)/bmass1
  return(bmass_change)
}

pdf("Change in biomass.pdf",width = 10,height = 3.5)
par(mfrow=c(1,3))
hist(sim_change(N_mean = 5),breaks=seq(-1,1,by=0.1), col="grey50", xlab="Proportional biomass change",main="Mean initial richness = 5")
abline(v=0,lty=2)
hist(sim_change(N_mean = 10),breaks=seq(-1,1,by=0.1),col="grey50", xlab="Proportional biomass change",main="Mean initial richness = 10")
abline(v=0,lty=2)
hist(sim_change(N_mean = 50),breaks=seq(-1,1,by=0.1),col="grey50", xlab="Proportional biomass change",main="Mean initial richness = 50")
abline(v=0,lty=2)
dev.off()

pdf("BEF Curve.pdf")
par(mfrow=c(1,1))
richV<-seq(0,55, length=1000)
plot(richV,4.33+richV^0.26, type='l', xlab="Species richness", ylab="Biomass", lwd=2)
abline(v=c(5,10,50), lty=2)
dev.off()

log(40/50)
log(4/5)
log(8/10
    )