#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure code
### Author: Mary O'Connor
### contributors: Jarrett Byrnes, Patrick Thompson
#######################################################################


library(lme4)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggExtra)

data <- read.csv("./data/SST5.csv")
data <- data[,-1]
data$TG1 <- factor(data$TG1)
mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

data$mod4.fitted <- fitted(mod4)

### FIGURE 1: 
### Extract slope coefficients for Figure 1B
## create the Entry, Study and Expt columns for later matching with random effect values
rand.cat <- ddply(dat, .(Entry, Study, ExptA, Sys1, TG1), summarize, mean(logY.rs))
names(rand.cat) <- c('Entry', 'Study', 'ExptA','System','TG1', 'meanlogY')
Entry.coefs <- data.frame(coef(mod4)$Entry) #this gives us the logS for entry, with study and expt level random effects (?)
Entry.coefs$Entry <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by = 'Entry', all = FALSE)

## test
ExptA.coefs <- data.frame(coef(mod4)$ExptA)
ExptA.coefs$ExptA <- rownames(ExptA.coefs)
Study.coefs <- data.frame(coef(mod4)$Study)
Study.coefs$Study <- rownames(Study.coefs)

S <- merge(S, ExptA.coefs[,c(2,12)], by = 'ExptA', all = FALSE)
S <- merge(S, Study.coefs[,c(2,12)], by = 'Study', all = FALSE)

##

b <- as.numeric(fixef(mod4)[2])
se <- function(x) sd(x)/sqrt(length(x))

## constructing predicted b values (slopes) from fixed effects for each group.
S$Trt <- paste(S$System, S$TG1)
S$Trt.term <- ifelse(S$Trt == 'A 1', (S$logSc.x + S$logSc.y + S$logSc - 2*b), 0)
S$Trt.term <- ifelse(S$Trt == 'A 2', (S$logSc.x + S$logSc.y + S$logSc - 2*b + S$logSc.TG12), S$Trt.term)
S$Trt.term <- ifelse(S$Trt == 'A 4', (S$logSc.x + S$logSc.y + S$logSc - 2*b + S$logSc.TG14), S$Trt.term)
S$Trt.term <- ifelse(S$Trt == 'T 1', (S$logSc.x + S$logSc.y + S$logSc - 2*b + S$logSc.Sys1T), S$Trt.term)
S$Trt.term <- ifelse(S$Trt == 'T 4', (S$logSc.x + S$logSc.y + S$logSc - 2*b + S$logSc.TG14 + S$logSc.Sys1T.TG14), S$Trt.term)


## producing Treatment mean coefficients
b.sum <- ddply(S, .(Trt), summarize, mean(Trt.term))
b.sum2 <- ddply(S, .(Trt), summarize, se(Trt.term))
b.sums <- merge(b.sum, b.sum2, by = 'Trt')
names(b.sums) <- c('group', 'est', 'se')

b.sums$group <- c("Aq. Primary Prod.", "Aq. Herbivores", "Aq. Detritovores", "Terr. Prim. Prod", "Terr. Detrit.")




### Figure 1A
fitted_plot <- ggplot() +
  geom_line(data=data, aes(x=exp(logSc+log(8)), y=exp(mod4.fitted), group=Entry), color="black", alpha = "0.3", size = .8) +
  theme_bw(base_size=17) +
  xlab("Species Richness (S)") +
  ylab("Estimated Biomass (Y)") +
  geom_text() +
  annotate("text", label = "A", x = 41, y = 15300, size = 6) + 
  removeGrid()

pdf(file = "Figure 1A.pdf", width = 4, height = 4)
fitted_plot 
dev.off()

### Figure 1B
pdf(file = "figure 1B.pdf", width = 4, height = 4)
par(
  oma = c(0,0,0,0), 
  mar = c(5,8,4,1), 
  mfcol = c(1,1)
)
 
plot(NULL,                                
     xlim = c(-0.2, 1.2),                        	
     ylim = c(.7, length(b.sums[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests.B <- (b.sums[,2])
ses.B <- (b.sums[,3])
var.names <- (b.sums[,1])

a <- 0
abline(v = 0, lty = 1, col = 'grey60')
for (i in 1:length(b.sums[,1])) {                                            
  lines(c(ests.B[i] + 1.96*ses.B[i], ests.B[i] - 1.96*ses.B[i]), c(i+a, i+a), col = 'gray60', lwd = 3)
  points(ests.B[i], i+a, pch = 19, cex = 1.5, col = 1) 
  text(-0.3, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names
  text(1.05, length(b.sums[,1])+ .1, 'B', cex = 1.8)
}

# add axes and labels
axis(side = 1)                                                                              
mtext(side = 1, "Slope coefficients", line = 3)                                             
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          

dev.off()


## Figure 1C
S.1 <- seq(1, 100, 1) #species richness before loss
S.2 = function(Y, S.1) exp((1/b) * ( log(Y) + b*log(S.1) ))
b <- 0.26

pdf(file = "figure 1C.pdf", width = 4, height = 4)
par(
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  mar = c(5,8,4,1),  # Inner margins are set through a little trial and error.
  mfcol = c(1,1)
)

plot(5*S.1^b ~ S.1, pch = '', xlim = c(1, 40), ylim = c(0, 80), axes = FALSE, xlab = 'Species richness (S)', ylab = 'Estimated Biomass (Y)')
axis(1, pos = 0)
axis(2, pos = 0, las = 2)
b <- 0.263
lines(5*S.1^b ~ S.1, lwd = 2, col = 1)

b <- 0.47
lines(5*S.1^b ~ S.1, lwd = 2, col = 'gray40')

b <- 0.52
lines(5*S.1^b ~ S.1, lwd = 2, col = 'gray60')

b <- -0.001
lines(5*S.1^b ~ S.1, lwd = 2, col = 'gray80')

legend(0, 80, c('Plants and Algae', 'Aq. Herbivores', 'Aq. Detritivores' ,'Terr. Detritivores'), pch=19, col = c('black', 'gray40', 'gray60', 'gray80'), bty = 'n')
text(37, 70, 'C', cex = 1.8)

dev.off()



### Figure 2
# simple plot
library(sjPlot)
library(arm)

pdf(file = "catplotEntry.pdf", width = 7.5, height = 3)
# make the caterpillar plot for Entry level coefs
cat.plot <- sjp.lmer(mod, sort = "logSc", ri.nr = 1, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

## percent of observations deviating from 0: 
test <- as.data.frame(cat.plot$plot[1])
test$coefnot0 <- ifelse(test$data.fade == FALSE, "0", "1")
length(test[(test$coefnot0 == "0"),9])/length(test$data.fade)

pdf(file = "catplotExptA.pdf", width = 7.5, height = 3)
# make the caterpillar plot for ExptA level coefs
cat.plot2 <- sjp.lmer(mod, sort = "logSc", ri.nr = 2, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

pdf(file = "catplotStudy.pdf", width = 7.5, height = 3)
# make the caterpillar plot for Study level coefs
cat.plot3 <- sjp.lmer(mod, sort = "logSc", ri.nr = 3, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()




### Figure 3


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
  xlim(0,2)+
  ylim(0,1.5)+
  xlab("Proportion of initial richness")+
  ylab("Proportion of initial function")+
  geom_vline(aes(xintercept = mean(RR,na.rm=T)),color="blue")+
  geom_hline(aes(yintercept = mean(RR,na.rm=T)),color="blue")+
  geom_hline(aes(yintercept = mean(YR,na.rm=T)),color="red")+
  geom_hline(aes(yintercept = 0.9), linetype = 2)+
  geom_vline(aes(xintercept = 0.65), linetype = 2)


ggMarginal(gg1,type = "histogram")        

geom_hline(aes(xintercept = 0.9,color=1, lin)+
             
             pdf("Simulated BEF change.pdf")
           ggMarginal(gg1,type = "histogram")
           dev.off()
           
           gg2<-ggplot(data.frame(RR=RR1, YR=RR1^b[2]),aes(x=RR,y=YR))+
             geom_point()+
             theme_bw(base_size = 16)+
             xlim(0,2)+
             ylim(0,1.5)+
             xlab("Proportion of initial richness")+
             ylab("Proportion of initial function")+
             geom_vline(aes(xintercept = mean(RR,na.rm=T)),color="blue")+
             geom_hline(aes(yintercept = mean(RR,na.rm=T)),color="blue")+
             geom_hline(aes(yintercept = mean(YR,na.rm=T)),color="red")+
             geom_hline(aes(yintercept = 0.9), linetype = 2)+
             geom_vline(aes(xintercept = 0.80), linetype = 2)
           
           ggMarginal(gg2, type = "histogram")
           
           
           
           
           