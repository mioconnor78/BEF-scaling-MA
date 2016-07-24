#######################################################################
### Toward predictive scaling coefficients for biodiversity and ecosystem functioning relationships

### O'Connor, Gonzalez, Byrnes, Cardinale, Duffy, Gamfeldt, Griffin, Hooper, Hungate, Paquette, Thompson, Dee and Matulich

### start with data processing march 2016 code

### Author of code: Mary O'Connor

#######################################################################

library(lme4)
library(MuMIn)
#library(AICcmodavg) # add this when you do m.avg
library(lmerTest)
#library(RLRsim)

# data <- read.csv("/Users/maryo/Documents/projects/BEF synthesis/BEF-scaling-MA/data/SST5.csv") 
data <- SST5

## equations for estimating degrees of freedom for models with different random effects, following Bolker et al wiki (ADD URL):
#q <- 2*2
#K <- function(x) length(fixef(x)) + (q*(q+1)/2)
#AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY.rs)/(length(data$logY.rs)-K(x)-1))


#########################################################################################
### 2. Hierarchical mixed effects model
########################################################################################

###### Part 1: The set of candiate level-1 models ##########
############################################################

## Determine best random effects structure for competing level-1 models (following O'Connor et al 2007)
## [use reml = FALSE for comparison with models with no random effects. Then, switch to reml = TRUE for random effects comparisons.]
## I'm not comparing lm with lmer fits due to differences in their estimation http://glmm.wikidot.com/faq

# Model 1 [Eqn 1a, main text], with different random effects structures
mod1 <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) +  (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod1i <- lmer(logY.rs ~ logSc + (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod1iii <- lmer(logY.rs ~ logSc + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod1v <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod1iv <- lmer(logY.rs ~ logSc + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)



# Model 2 [Eqn 2, main text]
mod2 <- lmer(logY.rs ~ logSc + log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod2i <- lmer(logY.rs ~ logSc + log(Tscale) +  (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod2iii <- lmer(logY.rs ~ logSc + log(Tscale) +  (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod2iv <- lmer(logY.rs ~ logSc + log(Tscale) +  (1 + logSc|Study) , data=data, REML = FALSE, na.action=na.omit)

mod2v <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)



# Model 3 [Eqn 3, main text]
mod3 <- lmer(logY.rs ~ logSc*log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod3i <- lmer(logY.rs ~ logSc*log(Tscale) +  (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod3iii <- lmer(logY.rs ~ logSc*log(Tscale) +  (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod3iv <- lmer(logY.rs ~ logSc*log(Tscale) +  (1 + logSc|Study) , data=data, REML = FALSE, na.action=na.omit)

mod3v <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

### Model comparison [Table A1]
## or using model.sel, without adjusting for degrees of freedom. appears that lmer now counts degrees of freedom as recommended by http://glmm.wikidot.com/faq. In previous explorations of these models, comparisons using adjusted likelihood ratio tests gave qualitatively the same outcomes. 
model.sel(mod1, mod1i, mod1iii, mod1iv, mod1v)
model.sel(mod2, mod2i, mod2iii, mod2iv, mod2v)
model.sel(mod3, mod3i, mod3iii, mod3iv, mod3v)


### [Table 2, MAIN TEXT] comparison of best level 1 model, using REML = FALSE 
mod1F <- lmer(logY.rs ~ logSc + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
mod2F <- lmer(logY.rs ~ logSc + log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
mod3F <- lmer(logY.rs ~ logSc*log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

model.sel(mod1F, mod2F, mod3F)

anova(mod2F, mod3F)
anova(mod3F, mod1F)

######################################################################
### Section 2: Testing different level-2 models. 
######################################################################

# Proceed here with random effects structure and level-1 model identifed above (Model 1 with full random effects).

# Model 4: Biological fixed factors that have been shown to not matter (system, trophic level) 
mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod4.2 <- lmer(logY.rs ~ logSc + log(Tscale) + logSc*Sys1 + logSc*TG1 + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Model 5: All biological fixed factors: ecosystem, trophic group, resource addition/reduction (adding Sys, TG, and res to the level-2 model)
mod5 <- lmer(logY.rs ~ logSc*Sys1*TG1  + logSc*restrt + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Model 6: All biological factors plus in interaction between richness and duration
mod6 <- lmer(logY.rs ~ logSc*Sys1*TG1 + logSc*restrt + logSc*log(maxTscale+1) + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Model 7: Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
mod7 <- lmer(logY.rs ~ logSc*restrt + logSc*log(maxTscale+1) + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Model 8: Experimental design factors (units, smax, time scale, Smax and units) 
mod8 <- lmer(logY.rs ~ logSc + log(Tscale) + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*log(maxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

# Model 9: Full model
mod9.1 <- lmer(logY.rs ~ logSc + log(Tscale) + logSc*Sys1 + logSc*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(maxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod9 <- lmer(logY.rs ~ log(Tscale) + logSc*Sys1*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(maxTscale+1) + (1 + logSc|Entry) + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)


###### Comparing models ##############
######################################

# [Table 3 Main Text] Model selection
model.sel(mod4, mod4.2, mod5, mod6, mod7, mod8, mod9.1, mod9, mod3F, mod2F)

# [Table 3 Main Text] Likelihood ratio tests
anova(mod4, mod5) # model 4 wins, so do next comparison
anova(mod4, mod6) # models not different, so proceed with model 4.2
anova(mod4, mod4.2) # models not different, so proceed with model 4.2
anova(mod4, mod9.2) 
anova(mod4, mod2F) 
anova(mod4, mod9) # can't do this one; models aren't nested. Will stop here; we're far enough down the list that I'm not interested in these comparisons anymore.
anova(mod4, mod3F)


# [Table A3: model summary] #use REML now for model coefficients
mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

summary(mod4)


## TABLE A2: does best model no longer need random effects? (needs them!)
mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod4i <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1|Entry) + (1|ExptA) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

mod4iii <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod4v <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod4iv <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)


model.sel(mod4, mod4i, mod4iii, mod4iv, mod4v)



##### 
## creating standardized residual plots for mixed models
#####
pdf(file = "Figure A2a.pdf", width = 4, height = 4)
plot(mod4, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 4")
dev.off()

pdf(file = "Figure A2b.pdf", width = 4, height = 4)
plot(mod4.2, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 4.2")
dev.off()

pdf(file = "Figure A2c.pdf", width = 4, height = 4)
plot(mod5, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 5")
dev.off()

pdf(file = "Figure A2d.pdf", width = 4, height = 4)
plot(mod6, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 6")
dev.off()

pdf(file = "Figure A2e.pdf", width = 4, height = 4)
plot(mod7, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 7")
dev.off()

pdf(file = "Figure A2f.pdf", width = 4, height = 4)
plot(mod8, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 8")
dev.off()

pdf(file = "Figure A2g.pdf", width = 4, height = 4)
plot(mod9, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 9")
dev.off()

pdf(file = "Figure A2h.pdf", width = 4, height = 4)
plot(mod9.1, ylab = "Standardized Residuals", xlab = "Fitted values", main = "Model 9.1")
dev.off()





######################################################################
### Section 3: Testing different level-2 models for final time point only, to test effect of duration without each year included. Same model set as above with Tscale term removed.
######################################################################

data <- SST5[(SST5$FinalT.x =='Y'),]

# Proceed here with random effects structure and level-1 model identifed above (Model 1 with full random effects).

# Model 4: Biological fixed factors that have been shown to not matter (system, trophic level) 
mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1  + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry) 

mod4.2 <- lmer(logY.rs ~ logSc + logSc*Sys1 + logSc*TG1 + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry) 

# Model 5: All biological fixed factors: ecosystem, trophic group, resource addition/reduction (adding Sys, TG, and res to the level-2 model)
mod5 <- lmer(logY.rs ~ logSc*Sys1*TG1  + logSc*restrt + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry) 

# Model 6: All biological factors plus in interaction between richness and duration
mod6 <- lmer(logY.rs ~ logSc*Sys1*TG1 + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry)

# Model 7: Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
mod7 <- lmer(logY.rs ~ logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry)

# Model 8: Experimental design factors (units, smax, time scale, Smax and units) 
mod8 <- lmer(logY.rs ~ logSc + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*log(MaxTscale+1) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) # + (1 + logSc|Entry)

# Model 9: Full model
mod9.1 <- lmer(logY.rs ~ logSc + logSc*Sys1 + logSc*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1) + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry) 

mod9 <- lmer(logY.rs ~ logSc*Sys1*TG1 + logSc*unit.types2 + logSc*Des1 + logSc*log(Smax) + logSc*restrt + logSc*log(MaxTscale+1)  + (1 + logSc|ExptA)  + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit) #+ (1 + logSc|Entry)

mod2 <- lmer(logY.rs ~ logSc  +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)


###### Comparing models ##############
######################################

# [Table S3] Model selection
model.sel(mod4, mod4.2, mod5, mod6, mod7, mod8, mod9.1, mod9, mod2)

# [Table 3 Main Text] Likelihood ratio tests
anova(mod4, mod5) # models not different, so proceed with model 4
anova(mod4, mod4.2) # models different, so proceed with model 4.2
anova(mod4.2, mod6) # models different, so proceed with model 6
anova(mod6, mod9) # models not different, so proceed with model 6
anova(mod6, mod2) # models different, so proceed with model 2
anova(mod2, mod8) # can't do this one; models aren't nested. Will stop here; we're far enough down the list that I'm not interested in these comparisons anymore.
anova(mod4, mod3F)

###################################
##### An alternate test for the effect of time within some groups (e.g., trophic groups)
### Table A4 ####

data <- read.csv("/Users/maryo/Documents/projects/BEF synthesis/BEF-scaling-MA/data/SST5.csv")

mod3F <- lmer(logY.rs ~ logSc*log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod3Fa <- lmer(logY.rs ~ logSc*log(Tscale)*Sys1 + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

mod3Fb <- lmer(logY.rs ~ logSc*log(Tscale)*TG1 + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

model.sel(mod3F, mod3Fa, mod3Fb)

summary(mod3Fb)
confint(mod3Fb)


##### exploring correlations between Smax and Smin
names(data)
plot(data$Smax ~ data$Slevels, ylab = 'Smax', xlab = 'Slevels', xlim = c(1, 10))
summary(lm(data$Smax ~ data$Slevels))

## create a column to identify minimum level of richness per entry
library(plyr)
tab1 <- ddply(data, .(Entry), summarize, min(richness))

## create a column to identify whether there are 2 or >2 levels of richness tested
data$no.levs <- ifelse(data$Slevels > 2, '2+', '2')
data$Entry <- as.numeric(as.character(data$Entry))

head(data)
data1 <- data[,-(3:46)]

# Create S from Figure 3 code
S$Entry <- as.numeric(as.character(S$Entry))
#data1$Entry <- as.numeric(as.character(data1$V1))
S1 <- merge(S, data1, by = 'Entry', all=FALSE, sort = TRUE)

tab2 <- ddply(S1, .(Entry, no.levs), summarise, min(slope))

pdf(file = "figure Slevels.pdf", width = 6, height = 4)
plot((tab2$..1)~as.factor(tab2$no.levs), xlab = 'No. of Richness levels tested in each Entry', ylab = 'Slope estimated by Model 4', )
dev.off()

dim(tab2[(tab2$no.levs == '2'),])


mean(S[(S$TG1=='1' & S$Syst == 'A'),]$slope)
mean(S[(S$TG1=='2'),]$slope)
mean(S[(S$TG1=='4'),]$slope)
mean(S[(S$Syst=='A'),]$slope)
