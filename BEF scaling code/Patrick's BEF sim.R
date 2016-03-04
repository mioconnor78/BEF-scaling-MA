######################################################################
#' @title Simulations of change in diversity leading to changes in function - alternate version
#' 
#' @author Patrick Thompson
#' @author patrickthompson9@gmail.com
#' 
#' @log
#'  2/20/2016 - First draft
######################################################################
require(dplyr)
require(data.table)
require(ggplot2)
require(tidyr)
require(ggExtra)

com<-5000 #number of communities to simulate
b <- c(0.25, 0.47, 0.53) # vector of scaling coefficients - not sure where the upper two come from but I assume that they correspond with herbivores and detritivores

RR1<-rnorm(com,mean=1,sd=0.4)
hist(RR1)
range(RR1)

RR1[RR1<0|RR1>2]<-NA
hist(RR1)

par(mfrow=c(1,4))
hist(RR1)# richness change histogram
hist(RR1^b[1], xlim=c(0,1.5)) #producer function
hist(RR1^b[2], xlim=c(0,1.5)) #herbivore function
hist(RR1^b[3], xlim=c(0,1.5)) #detritivore function

gg1<-ggplot(data.frame(RR=RR1,YR=RR1^b[1]),aes(x=RR,y=YR))+
  geom_point()+
  theme_bw(base_size = 16)+
  xlab("Proportion of initial richness")+
  ylab("Proportion of initial function")+
  geom_vline(aes(xintercept = mean(RR,na.rm=T)),color="blue")+
  geom_hline(aes(yintercept = mean(RR,na.rm=T)),color="blue")+
  geom_hline(aes(yintercept = mean(YR,na.rm=T)),color="red")

pdf("Simulated BEF change.pdf")
ggMarginal(gg1,type = "histogram")
dev.off()








#older sims####
RR1<-rnorm(com,mean=1,sd=0.1)
RR0.75<-rnorm(com,mean=0.75,sd=0.1)
RR0.5<-rnorm(com,mean=0.5,sd=0.1)
RR0.25<-rnorm(com,mean=0.25,sd=0.1)


BEF_change<-data.table(RR_mean=paste("S2/S1 =",c(rep(1,com),rep(0.75,com),rep(0.5,com),rep(0.25,com))),Richness=c(RR1,RR0.75,RR0.5,RR0.25))
BEF_change$Richness[BEF_change$Richness<0]<-0

BEF_change<-BEF_change%>%
  mutate(Producers=Richness^bV[1],Herbivores=Richness^bV[2],Carnivores=Richness^bV[3])

BEF_change<-gather(BEF_change,key = Type,value = Change,Richness:Carnivores)

BEF_change<-BEF_change%>%
  group_by(RR_mean,Type)%>%
  mutate(Mean_change=mean(Change))

pdf("Simulated BEF change by trophic level.pdf", width = 12,height=12)
ggplot(BEF_change,aes(x=Change))+
  geom_histogram(binwidth = 0.01)+
  facet_grid(Type~RR_mean,scale="free_y")+
  theme_bw(base_size = 16)+
  geom_vline(aes(xintercept=Mean_change), color="red", linetype=2)+
  geom_vline(xintercept = 1,col="blue",linetype=2)+
  xlab("Proportion of initial")
dev.off()


#even older sims###
b<-0.26

RR1<-rnorm(1000,mean=1,sd=0.1)
RR0.8<-rnorm(1000,mean=0.8,sd=0.1)
RR0.6<-rnorm(1000,mean=0.6,sd=0.1)

pdf("BEF scaling sim.pdf",width = 11,height = 8.5)
par(mfrow=c(2,3))
hist(RR1, xlab="Proportional richness change (S.2/S.1)", main="",breaks=seq(min(c(RR1,RR0.8,RR0.6)),max(c(RR1,RR0.8,RR0.6))+0.05,by=0.05))
abline(v=1, col=4, lwd=2)
abline(v=mean(RR1), col=2, lwd=2)

hist(RR0.8,xlab="Proportional richness change (S.2/S.1)", main="",breaks=seq(min(c(RR1,RR0.8,RR0.6)),max(c(RR1,RR0.8,RR0.6))+0.05,by=0.05))
abline(v=1, col=4, lwd=2)
abline(v=mean(RR0.8), col=2, lwd=2)

hist(RR0.6,xlab="Proportional richness change (S.2/S.1)", main="",breaks=seq(min(c(RR1,RR0.8,RR0.6)),max(c(RR1,RR0.8,RR0.6))+0.05,by=0.05))
abline(v=1, col=4, lwd=2)
abline(v=mean(RR0.6), col=2, lwd=2)

hist(RR1^b,xlab="Expected biomass change (Y.2/Y.1)", main="",breaks=seq(min(c(RR1^b,RR0.8^b,RR0.6^b)),max(c(RR1^b,RR0.8^b,RR0.6^b))+0.01,by=0.01))
abline(v=1, col=4, lwd=2)
abline(v=mean(RR1^b), col=2, lwd=2)

hist(RR0.8^b,xlab="Expected biomass change (Y.2/Y.1)", main="",breaks=seq(min(c(RR1^b,RR0.8^b,RR0.6^b)),max(c(RR1^b,RR0.8^b,RR0.6^b))+0.01,by=0.01))
abline(v=1, col=4, lwd=2)
abline(v=mean(RR0.8^b), col=2, lwd=2)

hist(RR0.6^b,xlab="Expected biomass change (Y.2/Y.1)", main="",breaks=seq(min(c(RR1^b,RR0.8^b,RR0.6^b)),max(c(RR1^b,RR0.8^b,RR0.6^b))+0.01,by=0.01))
abline(v=1, col=4, lwd=2)
abline(v=mean(RR0.6^b), col=2, lwd=2)
dev.off()

pdf("BEF scaling curve.pdf")
RR<-seq(0,2,length=1000)
par(mfrow=c(1,1))
plot(RR,(RR^b), type='l', xlab="S.2/S.1", ylab="Y.2/Y.1")
abline(v=range(RR1), col=2,lty=2)
abline(v=range(RR1), col=2,lty=2)
abline(v=range(RR0.8),col=3,lty=2)
abline(v=range(RR0.6),col=4,lty=2)
legend("bottomright",legend = c(1,0.8,0.6),col = 2:4,title="Proportion of\ninitial biomass", lty=2, bty='n')
dev.off()