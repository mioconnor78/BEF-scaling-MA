######################################################################
#' @title Simulations of change in diversity leading to changes in function - alternate version
#' 
#' @author Patrick Thompson
#' @author patrickthompson9@gmail.com
#' 
#' @log
#'  2/20/2016 - First draft
######################################################################

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

