#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### mixed effects analysis
### Jan 06 2015; Author: Mary O'Connor
#######################################################################

library(lme4)
library(MuMIn)
library(AICcmodavg)
library(RLRsim)
library(lmerTest)
#library(qpcR)

#########################################################################################
### 2. Hierarchical mixed effects model
#########################################################################################

###### The set of models ##########
###################################

data <- SST4

## determine best random effects structure for competing level-1 models (following O'Connor et al 2007)

# basic model (The level-1 model w/ time and logSc and time*logSc)  # add in fixed ints varying slopes
modBasic <- lmer(logYst ~ logSc*log(Tscale) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBasici <- lmer(logYst ~ logSc*log(Tscale) + (1|Entry) + (1|Study), data=data, REML = FALSE, na.action=na.omit)
modBasicii <- lm(logYst ~ logSc*log(Tscale), data=data, na.action=na.omit)

# basic2 (The level-1 model w/ time and logSc)
modBasic2 <- lmer(logYst ~ logSc+log(Tscale) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBasic2i <- lmer(logYst ~ logSc+log(Tscale) + (1|Entry) + (1|Study), data=data, REML = FALSE, na.action=na.omit)
modBasic2ii <- lm(logYst ~ logSc+log(Tscale), data=data, na.action=na.omit)

# basic3 (The level-1 model w/ logSc)
modBasic3 <- lmer(logYst ~ logSc + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBasic3i <- lmer(logYst ~ logSc + (1|Entry) + (1|Study), data=data, REML = FALSE, na.action=na.omit)
modBasic3ii <- lm(logYst ~ logSc, data=data, na.action=na.omit)



## compare with AIC adjusted for different degrees of freedome
q <- 2*2
K <- function(x) length(fixef(x)) + (q*(q+1)/2)
AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY)/(length(data$logY)-K(x)-1))
AICc.mem(modBasic) -> modBasic.AIC
AICc.mem(modBasic2) -> modBasic2.AIC
AICc.mem(modBasic3) -> modBasic3.AIC

#for fewer random effects
q <- 2 * 1 
AICc.mem(modBasici) -> modBasici.AIC
AICc.mem(modBasic2i) -> modBasic2i.AIC
AICc.mem(modBasic3i) -> modBasic3i.AIC

q <- 1 
AIC(modBasicii) -> modBasicii.AIC
AIC(modBasic2ii) -> modBasic2ii.AIC
AIC(modBasic3ii) -> modBasic3ii.AIC

## model comparison:
Basic <- cbind(modBasic.AIC, modBasici.AIC, modBasicii.AIC)
Basic2 <- cbind(modBasic2.AIC, modBasic2i.AIC, modBasic2ii.AIC)
Basic3 <- cbind(modBasic3.AIC, modBasic3ii.AIC, modBasic3ii.AIC)
Basic <- as.data.frame(rbind(Basic, cbind(logLik(modBasic), logLik(modBasici), logLik(modBasicii))))

#likelihood ratio test  ## HAVING SOME TROUBLE HERE
#1-pchisq(2*(logLik(mod1)-logLik(mod0)),1)
# D <- 2*(logLik(modBasic)-logLik(modBasici))
# D
# 1 - pchisq(D, 1)
# .5*(1-pchisq(LRstat,1)+1-pchisq(LRstat,2))

rand(modBasic)  # I wasn't believing the chi squared test b/c they were 0. but this shows that both ranefs are needed
#Analysis of Random effects Table:
 ##logSc:Entry     59      2   2e-13 ***
  #logSc:Study    115      2  <2e-16 ***

rand(modBasic2) 
#Analysis of Random effects Table:
#  Chi.sq Chi.DF p.value    
#logSc:Entry   64.9      2   8e-15 ***
#  logSc:Study  137.1      2  <2e-16 ***

rand(modBasic3) 

anova(modBasic, modBasic2)
anova(modBasic, modBasic3)

model.sel(modBasic, modBasic2, modBasic3)

### Set of candidate models #### 
# Full model 
modFM<-lmer(logYst ~ logSc:log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*unit.types2 + logSc*HigherT + logSc*logSmax + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# biological fixed factors that have been shown to not matter (system, trophic level, higher trophic level present) 
modBtrophic<-lmer(logYst ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
modBrt<-lmer(logYst ~ logSc*log(Tscale) + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# all biological fixed factors: ecosystem, trophic group, consumer presence, resource addition/reduction (adding Sys, TG, higherT, res and random effects to the level-2 model)
modBall<-lmer(logYst ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

modBallT<-lmer(logYst ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Experimental design factors (units, smax, time scale) [adding Duration.max, Smax and units to the level 2 model]
modExp<-lmer(logYst ~ logSc*log(Tscale) + logSc*unit.types2 + logSc*logSmax + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

modExp2<-lmer(logYst ~ logSc*log(Tscale) + logSc*log(Smax) + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

modExp3<-lmer(logYst ~ logSc*log(Tscale) + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

###### Comparing models ##############
######################################

model.sel(modFM, modBtrophic, modBrt, modBall, modExp, modBallT, modBasic) 


# In this section I am calculating AICc by hand, and accounting for degrees of freedom in the random effects according to Bolker et al. http://glmm.wikidot.com/faq
#AICc = -2*logLik(mod) + 2*K*(n/(n-K-1))
# the current code for model comparison below requires qpcR, which isn't working. I know from past analyses that the models with fewer random effects are terrible. So I'm going ahead with the model.avg command; it's imperfect but in this case the difference in df is not affecting the overall results.

q <- 2*2
K <- function(x) length(fixef(x)) + (q*(q+1)/2)
AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY)/(length(data$logY)-K(x)-1))
AICc.mem(modFM)
AIC.sum <- as.data.frame(cbind(AICc.mem(modFM), AICc.mem(modBtrophic), AICc.mem(modBrt), AICc.mem(modBall), AICc.mem(modExp), AICc.mem(modBasic)))
names(AIC.sum) <- c('modFM', 'modBtrophic', 'modBrt', 'modBall', 'modExp', 'modBasic')

#for fewer random effects
q <- 2 * 1
#AIC.sum$modBasici <- 
  AICc.mem(modBasici)

q <- 1
#AIC.sum$modBasicii <- 
  AIC(modBasicii)
#AIC.sum$modBasic3 <- AICc.mem(modBasic3)

AICc.mem(modBasic)

AIC.sum

## model averaging:
model.avg(modBtrophic, modBall, modBallT) -> m.avg  #modFM, 
m.avg

confint(m.avg)

summary(modBtrophic)

confint(modBtrophic) 


## save model objects for different datasets for Figure 2 plotting
data <- SST4
modBtrophic4<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBall<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBallT<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
m.avg.4 <- model.avg(modBtrophic4, modBall, modBallT)


data <- SST5
modBtrophic5<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBall<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBallT<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
m.avg.5 <- model.avg(modBtrophic5, modBall, modBallT)
