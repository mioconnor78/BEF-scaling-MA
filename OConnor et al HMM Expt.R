#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### mixed effects analysis, groups: ExptA and Study
### March 2015; Author: Mary O'Connor

### using this file for final analysis
#######################################################################

library(nlme)
library(MuMIn)
library(AICcmodavg)
library(RLRsim)
library(lmerTest)
library(bbmle)

#########################################################################################
### 2. Hierarchical mixed effects model
########################################################################################

###### The set of models ##########
###################################

data <- SST5
data$lnTscale <- log(data$Tscale)

## determine best random effects structure for competing level-1 models (following O'Connor et al 2007)

# basic model (The level-1 model w/ time and logSc and time*logSc) 
modBasic <- lme(logY.rs ~ logSc*lnTscale, random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML",  correlation = corCAR1(), na.action=na.omit, control=nlmeControl(minAbsParApVar=0.001, opt="nlminb", minScale=10e-10))      

modBasic1 <- lme(logY.rs ~ logSc*lnTscale, random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit, control=nlmeControl(minAbsParApVar=0.001, opt="nlminb", minScale=10e-10))                                                                                                                                                        

anova(modBasic, modBasic1)
AIC(modBasic, modBasic1)

modBasici <- lme(logY.rs ~ logSc*lnTscale, random = ~ 1 | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit) 
modBasicii <- lm(logY.rs ~ logSc*lnTscale, data=data, na.action=na.omit)

# basic2 (The level-1 model w/ time and logSc)
modBasic2 <- lme(logY.rs ~ logSc + lnTscale, random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit)
modBasic2i <- lme(logY.rs ~ logSc + lnTscale, random = ~ 1 | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit)
modBasic2ii <- lm(logY.rs ~ logSc+log(Tscale), data=data, na.action=na.omit)

# basic3 (The level-1 model w/ logSc)
modBasic3 <- lme(logY.rs ~ logSc, random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit)
modBasic3i <- lme(logY.rs ~ logSc, random = ~ 1 | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit)
modBasic3ii <- lm(logY.rs ~ logSc, data=data, na.action=na.omit)

# BASIC MODEL COMPARISON: 3 WAYS
#can't do likelihood ratio tests on the model with no variance components without switching to nlme, but lme4 is better for two random effects. using AIC and delta AIC is also accepted, and is convincing given the large differences between these models. 
# but, if I want to hang my hat on the presence of random effects, I could rewrite in nlme and to LR tests. 

## 1. Chi-square by hand: 
x2 = -2*logLik(modBasicii, REML=T) +2*logLik(modBasici, REML=T)
#pchisq(x2, df=1, lower.tail=F)
# for a comparison of a model with no random effects vs one with random intercept

0.5*(1-pchisq(x2, 1))

x2 = -2*logLik(modBasic, REML=T) +2*logLik(modBasici, REML=T)

# for a comparison of a model with no random effects vs one with random intercept
0.5*((1-pchisq(x2, 1)) + (1-pchisq(x2, 1)))


## 2. model comparison by AIC, without adusting degrees of freedom
bbmle::AICtab(modBasic, modBasici, modBasicii)
bbmle::AICtab(modBasic2, modBasic2i, modBasic2ii)
bbmle::AICtab(modBasic3, modBasic3i, modBasic3ii)


## 3. compare with AIC adjusted for different degrees of freedom
q <- 2*2
K <- function(x) length(fixef(x)) + (q*(q+1)/2)
AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY.rs)/(length(data$logY.rs)-K(x)-1))
AICc.mem(modBasic) -> modBasic.AIC
AICc.mem(modBasic2) -> modBasic2.AIC
AICc.mem(modBasic3) -> modBasic3.AIC

AICc.mem(modBtrophicii)
bbmle::AICtab(modBtrophic, modBtrophici, modBtrophicii)

#for fewer random effects
q <- 2 * 1 
AICc.mem(modBasici) -> modBasici.AIC
AICc.mem(modBasic2i) -> modBasic2i.AIC
AICc.mem(modBasic3i) -> modBasic3i.AIC

q <- 1 
AIC(modBasicii) -> modBasicii.AIC
AIC(modBasic2ii) -> modBasic2ii.AIC
AIC(modBasic3ii) -> modBasic3ii.AIC

Basic <- cbind(modBasic.AIC, modBasici.AIC, modBasicii.AIC)
Basic2 <- cbind(modBasic2.AIC, modBasic2i.AIC, modBasic2ii.AIC)
Basic3 <- cbind(modBasic3.AIC, modBasic3i.AIC, modBasic3ii.AIC)
#Basic.comb <- as.data.frame(rbind(Basic, cbind(logLik(modBasic), logLik(modBasici), logLik(modBasicii))))


rand(modBasic)  # I wasn't believing the chi squared test b/c they were 0. but this shows that both ranefs are needed

rand(modBasici) 
rand(modBasic2) 
rand(modBasic2i) 
rand(modBasic3) 
rand(modBasic3i) 

### having determined that random effects are needed in all models, compete basic models (Table 1)
anova(modBasic, modBasic2)
anova(modBasic2, modBasic3)

model.sel(modBasic, modBasic2, modBasic3)

# Full model 
modFull <- lme(logY.rs ~ logSc*lnTscale + logSc*Sys1  + logSc*TG1 + logSc*unit.types2 + logSc*HigherT + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1), random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML",   na.action=na.omit, control=nlmeControl(minAbsParApVar=0.001, opt="nlminb", minScale=10e-10))   

#correlation = corCAR1(),

modFull1 <- lme(logY.rs ~ logSc*lnTscale + logSc*Sys1  + logSc*TG1 + logSc*unit.types2 + logSc*HigherT + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1), random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML",  na.action=na.omit, control=nlmeControl(minAbsParApVar=0.001, opt="nlminb", minScale=10e-10))    

# biological fixed factors that have been shown to not matter (system, trophic level, higher trophic level present) 
modBtrophic<-lme(logY.rs ~ logSc*lnTscale + logSc*Sys1 + logSc*TG1 + logSc*HigherT, random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML", na.action=na.omit)

modBtrophici<-lme(logY.rs ~ logSc*lnTscale + logSc*Sys1 + logSc*TG1 + logSc*HigherT, random = ~ 1 | Study/ExptA/Entry, data=data, method = "REML", na.action=na.omit)

modBtrophicii<-lm(logY.rs ~ logSc*lnTscale + logSc*Sys1 + logSc*TG1 + logSc*HigherT, data=data, na.action=na.omit)

# Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
modBrt<-lme(logY.rs ~ logSc*lnTscale + logSc*restrt + logSc*log(MaxTscale+1), random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML", na.action=na.omit)

# all biological fixed factors: ecosystem, trophic group, consumer presence, resource addition/reduction (adding Sys, TG, higherT, res and random effects to the level-2 model)
modBall<-lme(logY.rs ~ logSc*lnTscale + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt, random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML", na.action=na.omit)

modBallT<-lme(logY.rs ~ logSc*lnTscale + logSc*Sys1 + logSc*TG1 + logSc*HigherT + logSc*restrt + logSc*log(MaxTscale+1), random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML", na.action=na.omit)

# Experimental design factors (units, smax, time scale) [adding Duration.max, Smax and units to the level 2 model]
modExp<-lme(logY.rs ~ logSc*lnTscale + logSc*unit.types2 + logSc*log(Smax) + logSc*log(MaxTscale+1), random = ~ logSc | Study/ExptA/Entry, data=data, method = "REML", na.action=na.omit)



###### Comparing models ##############
######################################




model.sel(modFull, modBtrophic, modBrt, modBall, modExp, modBasic1, modBallT)

## an attempt to get standardized coefficients 
## http://stackoverflow.com/questions/25142901/standardized-coefficients-for-lmer-model
stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

fixef(modBrt)
stdCoef.merMod(modBrt)

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
model.avg(modBtrophic, modBall) -> m.avg  #modFM, 
m.avg

confint(m.avg)

summary(modBtrophic)

confint(modBtrophic) 

### playing with caterpillar plots
resids <- residuals(modBtrophic)

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
