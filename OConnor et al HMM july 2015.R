#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### mixed effects analysis
### June 2015; Author: Mary O'Connor

### June 2015: FINAL ANALYSIS CODE, includes 3 random effects levels and standardized coefs, 3 way interaction
#######################################################################

library(lme4)
library(MuMIn)
library(AICcmodavg) # add this when you do m.avg
library(lmerTest)

data <- SST5

## equations for estimating degrees of freedom for models with different random effects, following Bolker et al wiki (ADD URL):
q <- 2*2
K <- function(x) length(fixef(x)) + (q*(q+1)/2)
AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY.rs)/(length(data$logY.rs)-K(x)-1))


#########################################################################################
### 2. Hierarchical mixed effects model
########################################################################################

###### Part 1: The set of candiate level-1 models ##########
############################################################

## Determine best random effects structure for competing level-1 models (following O'Connor et al 2007)
## [use reml = FALSE for comparison with models with no random effects. Then, switch to reml = TRUE for random effects comparisons.]

# Model 1 [Eqn 1, main text], with different random effects structures
mod1 <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod1i <- lmer(logY.rs ~ logSc +  (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod1iii <- lmer(logY.rs ~ logSc +  (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod1iv <- lmer(logY.rs ~ logSc  +  (1 + logSc|Study) , data=data, REML = FALSE, na.action=na.omit)

mod1ii <- lm(logY.rs ~ logSc, data=data, na.action=na.omit)

# Model 2 [Eqn 2, main text]
mod2 <- lmer(logY.rs ~ logSc + log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod2i <- lmer(logY.rs ~ logSc + log(Tscale) +  (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod2iii <- lmer(logY.rs ~ logSc + log(Tscale) +  (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod2iv <- lmer(logY.rs ~ logSc + log(Tscale) +  (1 + logSc|Study) , data=data, REML = FALSE, na.action=na.omit)

mod2ii <- lm(logY.rs ~ logSc + log(Tscale), data=data, na.action=na.omit)

# Model 3 [Eqn 3, main text]
mod3 <- lmer(logY.rs ~ logSc*log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod3i <- lmer(logY.rs ~ logSc*log(Tscale) +  (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod3ii <- lm(logY.rs ~ logSc*log(Tscale), data=data, na.action=na.omit)

mod3iii <- lmer(logY.rs ~ logSc*log(Tscale) +  (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod3iv <- lmer(logY.rs ~ logSc*log(Tscale) +  (1 + logSc|Study) , data=data, REML = FALSE, na.action=na.omit)


# BASIC MODEL COMPARISON: 3 WAYS
# can't do likelihood ratio tests on the model with no variance components without switching to nlme, but lme4 is better for two random effects. using AIC and delta AIC is also accepted, and is convincing given the large differences between these models. 
# but, if I want to hang my hat on the presence of random effects, I could rewrite in nlme and to LR tests. 

## 1. Chi-square by hand: 
x2 = -2*logLik(mod3ii, REML=T) +2*logLik(mod3i, REML=T)
#pchisq(x2, df=1, lower.tail=F)
# for a comparison of a model with no random effects vs one with random intercept

0.5*(1-pchisq(x2, 1))

x2 = -2*logLik(mod3, REML=T) +2*logLik(mod3i, REML=T)

# for a comparison of a model with no random effects vs one with random intercept
0.5*((1-pchisq(x2, 1)) + (1-pchisq(x2, 1)))


## 2. model comparison by AIC, without adusting degrees of freedom; bbmle is not available for R 3.2.1
#bbmle::AICtab(modBasic, modBasici, modBasicii, modBasiciii, modBasiciv)
#bbmle::AICtab(modBasic2, modBasic2i, modBasic2ii)
#bbmle::AICtab(modBasic3, modBasic3i, modBasic3ii)


## 3. compare with AIC adjusted for different degrees of freedom
q <- 3*2
K <- function(x) length(fixef(x)) + (q*(q+1)/2)
AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY.rs)/(length(data$logY.rs)-K(x)-1))
AICc.mem(mod3) -> mod3.AIC
AICc.mem(mod2) -> mod2.AIC
AICc.mem(mod1) -> mod1.AIC

#for fewer random effects
q <- 3 * 1 
AICc.mem(mod3i) -> mod3i.AIC
AICc.mem(mod2i) -> mod2i.AIC
AICc.mem(mod1i) -> mod1i.AIC

q <- 1 
AIC(mod3ii) -> mod3ii.AIC
AIC(mod2ii) -> mod2ii.AIC
AIC(mod1ii) -> mod1ii.AIC

q <- 3*1 
AIC(mod1iii) -> mod1iii.AIC
AIC(mod3iii) -> mod3iii.AIC
AIC(mod2iii) -> mod2iii.AIC

q <- 2*1
AIC(mod1iv) -> mod1iv.AIC
AIC(mod3iv) -> mod3iv.AIC
AIC(mod2iv) -> mod2iv.AIC

mod3.comp <- cbind(mod3.AIC, mod3i.AIC, mod3ii.AIC, mod3iii.AIC, mod3iv.AIC)
mod2.comp <- cbind(mod2.AIC, mod2i.AIC, mod2ii.AIC, mod2iii.AIC, mod2iv.AIC)
mod1.comp <- cbind(mod1.AIC, mod1i.AIC, mod1ii.AIC, mod1iii.AIC, mod1iv.AIC)

## or using model.sel, without adjusting for degrees of freedom [Table A1]
model.sel(mod1, mod1i, mod1ii, mod1iii, mod1iv)
model.sel(mod3, mod3i, mod3ii, mod3iii, mod3iv)
model.sel(mod2, mod2i, mod2ii, mod2iii, mod2iv)


### having determined that random effects are needed in all models, compete basic models (Table 1)
anova(modBasic, modBasic2)
anova(modBasic2, modBasic3)

### [Table 2, MAIN TEXT] comparison of best level 1 model, using REML = TRUE 
mod1F <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)
mod2F <- lmer(logY.rs ~ logSc + log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)
mod3F <- lmer(logY.rs ~ logSc*log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

model.sel(mod1F, mod2F, mod3F)

anova(mod1F, mod2F)
anova(mod1F, mod3F)

### Section 2: Testing different level 2 models. 
# biological fixed factors that have been shown to not matter (system, trophic level) 
# old name: modBtrophic
mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

# mod4.1 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod4.2 <- lmer(logY.rs ~ logSc + log(Tscale) + logSc*Sys1 + logSc*TG1 + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

#modB3trophic.so<-lmer(logY.rs ~ logSc*Sys1 + logSc*TG1 + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

# all biological fixed factors: ecosystem, trophic group, resource addition/reduction (adding Sys, TG, higherT, res and random effects to the level-2 model) (old names: modBall, modBallT)
mod5 <- lmer(logY.rs ~ logSc*Sys1*TG1  + logSc*restrt + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod6 <- lmer(logY.rs ~ logSc*Sys1*TG1 + logSc*restrt + logSc*log(MaxTscale+1) + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

# mod6.1 <- lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1*TG1 + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)


# Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
# old name: modBrt
# mod7.1 <- lmer(logY.rs ~ logSc*log(Tscale) + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod7 <- lmer(logY.rs ~ logSc*restrt + logSc*log(MaxTscale+1) + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

# Experimental design factors (units, smax, time scale) [adding Duration.max, Smax and units to the level 2 model]  #old name: modExp
mod8 <- lmer(logY.rs ~ logSc + log(Tscale) + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

# Full model
# old name: modFM
# mod9 <- lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod9.1 <- lmer(logY.rs ~ logSc + log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod9.2<-lmer(logY.rs ~ log(Tscale) + logSc*Sys1*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)



## does best model no longer need random effects? (needs them!)
modB3trophic2<-lmer(logY.rs ~ logSc+log(Tscale) + logSc*Sys1*TG1 + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

modB3trophic2a<-lmer(logY.rs ~ logSc+log(Tscale) + logSc*Sys1*TG1 + (1 + logSc|Entry), data=data, REML = FALSE, na.action=na.omit)

modB3trophic2b<-lmer(logY.rs ~ logSc+log(Tscale) + logSc*Sys1*TG1 + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

modBtrophicii<-lmer(logY.rs ~ logSc+log(Tscale) + logSc*Sys1*TG1 + (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

modBtrophiciii<-lm(logY.rs ~ logSc + logSc*Sys1*TG1, data=data, na.action=na.omit)

model.sel(modB3trophic2, modBtrophicii, modBtrophiciii, modB3trophic2a, modB3trophic2b)

###### Comparing models ##############
######################################

model.sel(mod4, mod4.2, mod5, mod6, mod7, mod8, mod9.1, mod9.2)


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
model.avg(modBallT, modBtrophic, modFM, modBall) -> m.avg  #modFM, 
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
