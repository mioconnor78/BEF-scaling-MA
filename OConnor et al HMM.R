#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### data formatting in preparation for mixed effects analysis
### Jan 06 2015; Author: Mary O'Connor
#######################################################################

library(lme4)
library(MuMIn)
library(AICcmodavg)
#library(qpcR)

#########################################################################################
### 2. Hierarchical mixed effects model
#########################################################################################

# for all models; 
# Levels: Study (i), Experiment (j) and Plot (k)
# Predictors (level): System (study), TG (expt), units (expt), higherT (expt), Smax (expt), logS (plot), Time (plot)

# Level 1 model (Eqn 1a from paper)
# ln(Y) = B.0ijk + B.1ijk*log(S.ijk) + B.2*log(Duration).ij + B.3*log(S.ijk).ijk*log(Duration).ijk + e.ijk

# Level 2 model (Eqn 1b from paper)
# B.0ij = y.00 + y.01*Sys.i + y.02*TG.ij + y.03*Units.ij + y.04*HT.ij + y.05*Smax.ij + y.06*Res.ij + y.07*max(lnT).ij + u.0i + u.0j
# B.1ij = y.10 + y.11*Sys.i + y.12*TG.ij + y.13*Units.ij + y.14*HT.ij + y.15*Smax.ij + y.16*Res.ij + y.17*max(lnT).ij + u.1i + u.1j

# Level 2 substituted back into Level 1 for full model: 
# ln(Y) = (y.00 + y.01*Sys.i + y.02*TG.ij + y.03*Units.ij + y.04*HT.ij + y.05*Smax.ij + y.06*Res.ij + y.07*max(lnT).ij + u.0i + u.0j).k 
#  + (y.10 + y.11*Sys.i + y.12*TG.ij + y.13*Units.ij + y.14*HT.ij + y.15*Smax.ij + y.16*Res.ij + y.17*max(lnT).ij + u.1i + u.1j)k*log(S.ijk) 
#  + B.2*log(Duration).ij + B.3*log(S.ijk).ijk*log(Duration).ijk + e.ijk

# Full model, rearranged:



###### The set of models ##########
###################################

data <- SST5

# Full model 
modFM<-lmer(logY.rs ~ logS*log(Tscale) + logS*Sys1  + logS*TG1 + logS*unit.types2 + logS*HigherT + logS*log(Smax) + logS*restrt + logS*log(MaxTscale+1) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)
modFMi<-lmer(logY ~ logS*log(Tscale) + logS*Sys1  + logS*TG1 + logS*unit.types2 + logS*HigherT + logS*log(Smax) + logS*restrt + logS*log(MaxTscale+1) + (1|Entry) + (1|Study), data=data, REML = FALSE, na.action=na.omit)
modFMii<-lm(logY ~ logS*log(Tscale) + logS*Sys1  + logS*TG1 + logS*unit.types2 + logS*HigherT + logS*log(Smax) + logS*restrt + logS*log(MaxTscale+1), data=data, na.action=na.omit)

# biological fixed factors that have been shown to not matter (system, trophic level, higher trophic level present) 
modBtrophic<-lmer(logY.rs ~ logS*log(Tscale) + logS*Sys1 + logS*TG1 + logS*HigherT + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
modBrt<-lmer(logY.rs ~ logS*log(Tscale) + logS*restrt + logS*log(MaxTscale+1) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# all biological fixed factors: ecosystem, trophic group, consumer presence, resource addition/reduction (adding Sys, TG, higherT, res and random effects to the level-2 model)
modBall<-lmer(logY.rs ~ logS*log(Tscale) + logS*Sys1 + logS*TG1 + logS*HigherT + logS*restrt + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# Experimental design factors (units, smax, time scale) [adding Duration.max, Smax and units to the level 2 model]
modExp<-lmer(logY.rs ~ logS*log(Tscale) + logS*unit.types2 + logS*log(Smax) + logS*log(MaxTscale+1) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# basic model (The level-1 model level-2 that includes with random effects only)
modBasic <- lmer(logY.rs ~ logS*log(Tscale) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

#modBasic1 <- lmer(logY ~ logS*log(Tscale) + (1|Mno) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

#modBasic2 <- lmer(logY ~ logS*log(Tscale) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

#modBasic3 <- lm(logY ~ logS*log(Tscale), data = data, na.action = na.omit)



###### Comparing models ##############
######################################

model.sel(modFM, modBtrophic, modBrt, modBall, modExp, modBasic, modFMi, modFMii)
# In this section I am calculating AICc by hand, and accounting for degrees of freedom in the random effects according to Bolker et al. http://glmm.wikidot.com/faq
#AICc = -2*logLik(mod) + 2*K*(n/(n-K-1))

q <- 2*2
K <- function(x) length(fixef(x)) + (q*(q+1)/2)
AICc.mem <- function(x) -2*as.numeric(logLik(x)) + 2*K(x)*(length(data$logY)/(length(data$logY)-K(x)-1))
AICc.mem(modFM)

AIC.sum <- as.data.frame(cbind(AICc.mem(modFM), AICc.mem(modBtrophic), AICc.mem(modBrt), AICc.mem(modBall), AICc.mem(modExp), AICc.mem(modBasic)))
names(AIC.sum) <- c('modFM', 'modBtrophic', 'modBrt', 'modBall', 'modExp', 'modBasic')

#for fewer random effects
q <- 2 * 1
AIC.sum$modFMi <- AICc.mem(modFMi)

q <- 1
AIC.sum$modFMii <- AIC(modFMii)
#AIC.sum$modBasic3 <- AICc.mem(modBasic3)

AIC.sum

summary(modBtrophic)
confint(modBtrophic)

## model averaging:
model.avg(modBtrophic, modBall) -> m.avg  #modFM, 
m.avg

confint(m.avg)

summary(modBtrophic)

confint(modBtrophic) 

### Figure 2 ####
estimates <- as.data.frame(m.avg[3])
var.names<-rownames(estimates)
var.names<-c('Intercept', 'logS', 'Ecosystem - Aq', 'TG: Herbivore', 'TG: Carnivore', 'TG: Detritovore', 'Consumer present', '+ Resource', '- Resource', 'log(Time in gen)', 'logS*Ecosystem', 'logS*Herbivore', 'logS*Carnivore', 'logS*Detritovore', 'logS*Consumer pres', 'logS * + Resource', 'logS* - Resource', 'logS*log(Time in gen)', 'units Density', 'units % Cover', 'log(Smax)', 'logS*Density', 'logS*mass', 'logS*log(Smax)')


## plotting from http://www.carlislerainey.com/Blog_Files/Blog_CoefficientPlots.R
# set the graphical parameters
par(
  family = "serif",  # I don't plot in anything but serif
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  mar = c(5,10,4,2)  # Inner margins are set through a little trial and error.
)

# create an empty plot for total customization
plot(NULL,                                	# create empty plot
     xlim = c(-3, 4),                      		# set xlim by guessing
     ylim = c(.7, length(estimates[,1]) + .3), 	# set ylim by the number of variables
     axes = F, xlab = NA, ylab = NA)       		# turn off axes and labels

# add the data
est <- estimates[,1] #[-1]                                             # conveniently store the estimates (minus the constant)
se <- estimates[,2]                                        		# conveniently store the std. errors (minus the constant)
for (i in 1:length(est)) {                                            # loop over a counter the length of the estimate vector
  points(est[i], i, pch = 19, cex = .75)                               # add the points to the plot
  #lines(c(est[i] + 1.64*se[i], est[i] - 1.64*se[i]), c(i, i))         # add the 90% confidence intervals (1.64)
  lines(c(estimates$avg.model.Upper.CI[i], estimates$avg.model.Lower.CI[i]), c(i, i))
  #lines(c(est[i] + .67*se[i], est[i] - .67*se[i]), c(i, i), lwd = 3)  # add the 50% confidence intervals
  text(-4.5, i, var.names[i], xpd = T, cex = .8)                      # add the variable names
}

# add axes and labels
axis(side = 1)                                                                                          # add bottom axis
abline(v = 0, lty = 3, col = "grey")                                                                    # add verticle line
mtext(side = 1, "Model-Averaged Coefficient", line = 3)                                              # label bottom axis
mtext(side = 3, "Model-averaged coefficients of\n logY = f(logS)", line = 1)   # add title
box()                                                                                                   # add lines around the plot

##############################################
#### plot just slope coefficients...
###############################################

estimates <- as.data.frame(m.avg[3])
var.names<-rownames(estimates)
## varnames SST4 w/ preds
var.names<-c('Intercept', 'ln(S)', 'ln(Duration)',  'Ecosystem - T', 'TG: Herbivore', 'TG: Predator', 'TG: Detritivore', '+ Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator', 'ln(S)*Detritivore', 'ln(S)* +Consumer',  '+Resource', '-Resource',  'ln(S) * +Resource', 'ln(S)* -Resource')
estimates$var.names<-var.names
sl.est <- estimates[-1,]
sl.est <- sl.est[-(3:7),]
sl.est <- sl.est[-(9:10),]


## varnames SST5 (w/o preds)
estimates <- as.data.frame(m.avg[3])
var.names<-c('Intercept', 'ln(S)', 'ln(Duration)',  'Ecosystem - Aq', 'TG: Herbivore', 'TG: Detritivore', '+ Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore', 'ln(S)* +Consumer',  '+Resource', '-Resource',  'ln(S) * +Resource', 'ln(S)* -Resource')
estimates$var.names<-var.names
sl.est <- estimates[-1,]
sl.est <- sl.est[-(3:6),]
sl.est <- sl.est[-(8:9),]
# add blank line for where preds would be
sl.est <- rbind(sl.est[1:5,], c('','','','','',''), sl.est[6:nrow(sl.est),])
row.names(sl.est) <- c('ln(S)', 'ln(Duration)',  'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', '','ln(S)*Detritivore', 'ln(S)* +Consumer', 'ln(S) * +Resource', 'ln(S)* -Resource')


## plotting from http://www.carlislerainey.com/Blog_Files/Blog_CoefficientPlots.R
# set the graphical parameters
par(
  family = "serif",  # I don't plot in anything but serif
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  mar = c(5,10,4,2)  # Inner margins are set through a little trial and error.
  #mfrow = c(1,1)
)

# create an empty plot for total customization
plot(NULL,                              		# create empty plot
     xlim = c(-0.6, 0.6),                      		# set xlim by guessing
     ylim = c(.7, length(sl.est[,1]) + .3), 	# set ylim by the number of variables
     axes = F, xlab = NA, ylab = NA)       		# turn off axes and labels

# add the data
est <- sl.est[,1]   # conveniently store the estimates (minus the constant)
abline(v = 0, lty = 3, lwd = 3, col = "grey")                                             
#se <- estimates[,2]      # conveniently store the std. errors (minus the constant)
for (i in 1:length(est)) {      # loop over a counter the length of the estimate vector
  points(est[i], i, pch = 19, cex = 1)           
  #lines(c(est[i] + 1.64*se[i], est[i] - 1.64*se[i]), c(i, i))   # add the 90% confidence intervals (1.64)
  lines(c(sl.est$avg.model.Upper.CI[i], sl.est$avg.model.Lower.CI[i]), c(i, i))
  text(-1, i, sl.est$var.names[i], xpd = T, cex = 0.9, adj = c(0.5,0.5))                      
}

# add axes and labels
axis(side = 1, cex.axis = 0.8)                                                                                          # add bottom axis
# abline(v = 0, lty = 3, lwd = 3, col = "grey")                                                                    # add verticle line
mtext(side = 1, "Model-Averaged Coefficients", line = 3, cex = 1.2)                                # label bottom axis
# mtext(side = 3, paste("Model-averaged Coefficients\n for terms in the BEFR", expression(beta[1][ij]), '+', expression(Beta[3][ij])), line = 1)   # add title
mtext(side = 3, "B. Predator Studies Removed", line = 1, cex = 1.2)   # add title
box()                                      

# A. Full dataset
# B. Predator Studies Removed



