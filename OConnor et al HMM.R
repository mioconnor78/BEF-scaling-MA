#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### data formatting in preparation for mixed effects analysis
### Jan 06 2015; Author: Mary O'Connor
#######################################################################

library(lme4)
library(MuMIn)
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

data <- SST2

# Full model 
modFM<-lmer(logY ~ logS*log(Tscale) + logS*Sys1  + logS*TG1 + logS*unit.types2 + logS*HigherT + logS*log(Smax) + logS*restrt + logS*log(MaxTscale+1) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)
modFMi<-lmer(logY ~ logS*log(Tscale) + logS*Sys1  + logS*TG1 + logS*unit.types2 + logS*HigherT + logS*log(Smax) + logS*restrt + logS*log(MaxTscale+1) + (1|Entry) + (1|Study), data=data, REML = FALSE, na.action=na.omit)
modFMii<-lm(logY ~ logS*log(Tscale) + logS*Sys1  + logS*TG1 + logS*unit.types2 + logS*HigherT + logS*log(Smax) + logS*restrt + logS*log(MaxTscale+1), data=data, na.action=na.omit)

# biological fixed factors that have been shown to not matter (system, trophic level, higher trophic level present) 
modBtrophic<-lmer(logY ~ logS*log(Tscale) + logS*Sys1 + logS*TG1 + logS*HigherT + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# Fixed factors that have been shown to matter (adding time, nutrients to level-2 model)
modBrt<-lmer(logY ~ logS*log(Tscale) + logS*restrt + logS*log(MaxTscale+1) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# all biological fixed factors: ecosystem, trophic group, consumer presence, resource addition/reduction (adding Sys, TG, higherT, res and random effects to the level-2 model)
modBall<-lmer(logY ~ logS*log(Tscale) + logS*Sys1 + logS*TG1 + logS*HigherT + logS*restrt + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# Experimental design factors (units, smax, time scale) [adding Duration.max, Smax and units to the level 2 model]
modExp<-lmer(logY ~ logS*log(Tscale) + logS*unit.types2 + logS*log(Smax) + logS*log(MaxTscale+1) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

# basic model (The level-1 model level-2 that includes with random effects only)
modBasic <- lmer(logY ~ logS*log(Tscale) + (1 + logS|Entry) + (1 + logS|Study), data=data, REML = FALSE, na.action=na.omit)

#modBasic1 <- lmer(logY ~ logS*log(Tscale) + (1|Mno) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

#modBasic2 <- lmer(logY ~ logS*log(Tscale) + (1|Study), data=data, REML = FALSE, na.action=na.omit)

#modBasic3 <- lm(logY ~ logS*log(Tscale), data = data, na.action = na.omit)



###### Comparing models ##############
######################################

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

# calculating weights and delta AIC
#round(Weights(AICc(modFM, modBtrophic, modBrt, modBall, modExp, modBasic)))
model.sel(modFM, modBtrophic, modBrt, modBall, modExp, modBasic)

#weights <- akaike.weights(AIC.sum) #requires qpcR, which isn't working
#results <- as.data.frame(rbind(weights$deltaAIC, weights$rel.LL, weights$weights))
#names(results) <- c('modFM', 'modBtrophic', 'modBrt', 'modBall', 'modExp', 'modBasic', 'modFMi') #, 'modFMii'
#rownames(results) <- c('deltaAIC', 'relLL','weights')
results

#SST2 (after code merge). different results...
            logLik    AICc   delta weight
modFM       -1114.846 2301.1  0.00 1     
modBtrophic -1137.202 2316.9 15.79 0     
modBall     -1135.115 2321.0 19.83 0     
modBasic    -1165.581 2353.3 52.17 0     
modExp      -1155.863 2354.3 53.11 0     
modBrt      -1163.294 2360.9 59.80 0       

#SST2
> results (before code merge)
modFM modBtrophic       modBrt   modBall       modExp     modBasic       modFMi  modFMii
deltaAIC 9.492287829   0.0000000 1.771517e+01 3.6930563 2.358546e+01 10.692942333 2.342200e+02 4342.987
relLL    0.008685121   1.0000000 1.422980e-04 0.1577840 7.559315e-06  0.004764936 1.379683e-51    0.000
weights  0.007414410   0.8536911 1.214785e-04 0.1346988 6.453320e-06  0.004067783 1.177823e-51    0.000


#SST2nc
modFM modBtrophic       modBrt   modBall       modExp     modBasic       modFMi  modFMii
deltaAIC 7.49927119   0.0000000 1.897615e+01 3.5123883 2.182002e+01 12.224708046 1.985566e+02 4315.648
relLL    0.02352632   1.0000000 7.574994e-05 0.1727009 1.827443e-05  0.002215330 7.655599e-44    0.000
weights  0.01962920   0.8343508 6.320203e-05 0.1440931 1.524728e-05  0.001848362 6.387456e-44    0.000






model.avg(modBtrophic, modBall) -> m.avg  #modFM, 

Call:
  model.avg.default(object = modBtrophic, modBall)

Component models: 
  ‘1/2/3/5/6/7/8/10/11’     ‘1/2/3/4/5/6/7/8/9/10/11’

Coefficients: 
  (Intercept)             logS      log(Tscale)             Sys1             TG12             TG13             TG14         HigherTY 
5.393812723      0.188949067     -0.005001719     -0.923645493     -1.535731706      0.317716990     -0.971956144     -0.235696176 
log(Tscale):logS        logS:Sys1        logS:TG12        logS:TG13        logS:TG14    HigherTY:logS       restrtincr        restrtred 
0.030255839      0.088800715      0.152902403     -0.303286499      0.066875787      0.021028595      0.186423211     -0.534005118 
logS:restrtincr   logS:restrtred 
-0.013008265     -0.075334776 
> confint(m.avg)
2.5 %      97.5 %
  (Intercept)       4.841334172  5.94629127
logS              0.101470097  0.27642804
log(Tscale)      -0.141891715  0.13188828
Sys1             -1.749689109 -0.09760188
TG12             -2.416489626 -0.65497379
TG13             -1.412528428  2.04796241
TG14             -2.449345785  0.50543350
HigherTY         -0.764572105  0.29317975
log(Tscale):logS  0.005903297  0.05460838
logS:Sys1        -0.041347946  0.21894937
logS:TG12        -0.029441297  0.33524610
logS:TG13        -0.595874969 -0.01069803
logS:TG14        -0.159874529  0.29362610
HigherTY:logS    -0.070576414  0.11263360
restrtincr       -0.060417097  0.43326352
restrtred        -1.563966609  0.49595637
logS:restrtincr  -0.065772764  0.03975623
logS:restrtred   -0.258997042  0.10832749



> confint(modBtrophic) #for nc dataset
Computing profile confidence intervals ...
2.5 %      97.5 %
  .sig01            0.865944654  0.99065108
.sig02           -0.336965667 -0.06703099
.sig03            0.104264253  0.14989400
.sig04            1.310848054  1.87638920
.sig05           -0.190321903  0.38353087
.sig06            0.138845702  0.22532603
.sigma            0.171615809  0.19219702
(Intercept)       4.859911524  5.99317578
logS              0.098037492  0.25906306
log(Tscale)      -0.139413186  0.15443089
Sys1             -1.908792384 -0.16557751
TG12             -2.425629034 -0.62367844
TG14             -2.464870169  0.57282053
HigherTY         -0.769443687  0.31361959
logS:log(Tscale) -0.005240245  0.04360586
logS:Sys1        -0.032189069  0.20986242
logS:TG12         0.007243528  0.35288929
logS:TG14        -0.115569499  0.29987000
logS:HigherTY    -0.052932534  0.11895928



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

#### could plot just slope coefficients...
###############################################

estimates <- as.data.frame(m.avg[3])
var.names<-rownames(estimates)
## varnames SST2
var.names<-c('Intercept', 'ln(S)', 'ln(Duration)',  'Ecosystem - Aq', 'TG: Herbivore', 'TG: Predator', 'TG: Detritivore', '+ Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator', 'ln(S)*Detritivore', 'ln(S)* + Consumer',  '+ Resource', '- Resource',  'ln(S) * + Resource', 'ln(S)* - Resource')
estimates$var.names<-var.names
sl.est <- estimates[-1,]
sl.est <- sl.est[-(3:7),]
sl.est <- sl.est[-(9:10),]


## varnames SST2nc
var.names<-c('Intercept', 'ln(S)', 'ln(Duration)',  'Ecosystem - Aq', 'TG: Herbivore', 'TG: Detritivore', '+ Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore', 'ln(S)* + Consumer',  '+ Resource', '- Resource',  'ln(S) * + Resource', 'ln(S)* - Resource')
estimates$var.names<-var.names
sl.est <- estimates[-1,]
sl.est <- sl.est[-(3:6),]
sl.est <- sl.est[-(8:9),]

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
est <- sl.est[,1] #[-1]   
abline(v = 0, lty = 3, lwd = 3, col = "grey")                                             # conveniently store the estimates (minus the constant)
#se <- estimates[,2]                                        		# conveniently store the std. errors (minus the constant)
for (i in 1:length(est)) {                                            # loop over a counter the length of the estimate vector
  points(est[i], i, pch = 19, cex = 1)                               # add the points to the plot
  #lines(c(est[i] + 1.64*se[i], est[i] - 1.64*se[i]), c(i, i))         # add the 90% confidence intervals (1.64)
  lines(c(sl.est$avg.model.Upper.CI[i], sl.est$avg.model.Lower.CI[i]), c(i, i))
  #lines(c(est[i] + .67*se[i], est[i] - .67*se[i]), c(i, i), lwd = 3)  # add the 50% confidence intervals
  text(-1, i, sl.est$var.names[i], xpd = T, cex = 1.2, adj = c(0.5,0.5))                      # add the variable names
}

# add axes and labels
axis(side = 1)                                                                                          # add bottom axis
#abline(v = 0, lty = 3, lwd = 3, col = "grey")                                                                    # add verticle line
mtext(side = 1, "Model-Averaged Coefficients", line = 3, cex = 1.2)                                # label bottom axis
#mtext(side = 3, paste("Model-averaged Coefficients\n for terms in the BEFR", expression(beta[1][ij]), '+', expression(Beta[3][ij])), line = 1)   # add title
mtext(side = 3, "B. Predator Studies Removed", line = 1, cex = 1.2)   # add title
box()                                      







#### summary of the best model: ModF1
> summary(modF1)
Linear mixed model fit by maximum likelihood ['lmerMod']
Formula: logY ~ logS * Sys1 + logS * TG1 + logS * HigherT + logS * log(Tscale) +      (1 + logS | Mno) + (1 + logS | Study) 
Data: data 

AIC       BIC    logLik  deviance 
2167.066  2282.155 -1062.533  2125.066 

Random effects:
  Groups   Name        Variance Std.Dev. Corr 
Mno      (Intercept) 0.84596  0.9198        
logS        0.01624  0.1275   -0.20
Study    (Intercept) 2.35819  1.5356        
logS        0.04377  0.2092   0.17 
Residual             0.03333  0.1826        
Number of obs: 1773, groups: Mno, 559; Study, 90

Fixed effects:
  Estimate Std. Error t value
(Intercept)       5.397708   0.281406  19.181
logS              0.188052   0.044462   4.230
Sys1             -0.925327   0.421490  -2.195
TG12             -1.535717   0.449546  -3.416
TG13              0.326702   0.882771   0.370
TG14             -0.967673   0.754008  -1.283
HigherTY         -0.237471   0.269809  -0.880
log(Tscale)      -0.004844   0.069866  -0.069
logS:Sys1         0.089479   0.066343   1.349
logS:TG12         0.152894   0.093031   1.643
logS:TG13        -0.303618   0.149250  -2.034
logS:TG14         0.066833   0.115680   0.578
logS:HigherTY     0.021382   0.046706   0.458
logS:log(Tscale)  0.030209   0.012424   2.432

Correlation of Fixed Effects:
  (Intr) logS   Sys1   TG12   TG13   TG14   HghrTY lg(Ts) lgS:S1 lS:TG12 lS:TG13 lS:TG14 lS:HTY
logS         0.026                                                                                       
Sys1        -0.522 -0.012                                                                                
TG12        -0.228 -0.006 -0.147                                                                         
TG13        -0.119 -0.008 -0.106  0.092                                                                  
TG14        -0.202  0.006 -0.135  0.158  0.107                                                           
HigherTY    -0.638  0.038  0.235  0.303  0.015  0.134                                                    
log(Tscale)  0.198 -0.035 -0.211 -0.151 -0.126 -0.188 -0.061                                             
logS:Sys1   -0.009 -0.545  0.054  0.022  0.004 -0.032 -0.027  0.040                                      
logS:TG12   -0.010 -0.138  0.021 -0.073 -0.011 -0.004  0.016  0.000 -0.249                               
logS:TG13   -0.008 -0.094  0.006 -0.014  0.036  0.002 -0.011 -0.003 -0.146  0.118                        
logS:TG14    0.007 -0.242 -0.032 -0.008  0.002  0.065 -0.014  0.035 -0.064  0.156   0.109                
logS:HghrTY  0.035 -0.701 -0.022  0.011 -0.010 -0.012 -0.050  0.028  0.299  0.157   0.038   0.163        
lgS:lg(Tsc) -0.030  0.250  0.031  0.007  0.000  0.028  0.026 -0.091 -0.285 -0.094  -0.097  -0.258  -0.104


#confint.merMod(modFM)  # returns the same values as confint(modFM)
Computing profile confidence intervals ...
2.5 %      97.5 %
  .sig01            0.861434181  0.98458813
.sig02           -0.327564241 -0.05828934
.sig03            0.104146030  0.15056943
.sig04            1.294743923  1.83681362
.sig05           -0.104616163  0.42480762
.sig06            0.168414312  0.25856906
.sigma            0.172559139  0.19347406
(Intercept)       4.844401981  5.95897336
logS              0.099001850  0.27565299
Sys1             -1.764606517 -0.09497573
TG12             -2.425821992 -0.65130768
TG13             -1.433824995  2.06635157
TG14             -2.465628734  0.52182158
HigherTY         -0.771287951  0.29250381
log(Tscale)      -0.144909871  0.13751441
logS:Sys1        -0.041086807  0.22221434
logS:TG12        -0.031107828  0.33621777
logS:TG13        -0.603120800 -0.01024197
logS:TG14        -0.160330592  0.29963184
logS:HigherTY    -0.070598814  0.11411579
logS:log(Tscale)  0.005027134  0.05484957



###### Calculating R2, based on nakagawa and schielzeth 
###### the sample code online includes the code mod@X, but this doesn't appear to work in the new lme4. I have figured out (!!) that the equiavalent that works is model.matrix(mod)

## January 21: I now notice that nakagawa adn schielzeth say this doesn't apply to random slopes. so i shoudl be including those components. They recommend estimating variance for the equivalent, random intercetps only model. 

mod<-modFMi
## for model F1
Fixed <- fixef(mod)[2]*model.matrix(mod)[,2] + fixef(mod)[3]*model.matrix(mod)[,3] + fixef(mod)[4]*model.matrix(mod)[,4] + fixef(mod)[5]*model.matrix(mod)[,5] + fixef(mod)[6]*model.matrix(mod)[,6] + fixef(mod)[7]*model.matrix(mod)[,7] + fixef(mod)[8]*model.matrix(mod)[,8] + fixef(mod)[9]*model.matrix(mod)[,9] + fixef(mod)[10]*model.matrix(mod)[,10] + fixef(mod)[11]*model.matrix(mod)[,11] + fixef(mod)[12]*model.matrix(mod)[,12] + fixef(mod)[13]*model.matrix(mod)[,13] + fixef(mod)[14]*model.matrix(mod)[,14] + fixef(mod)[15]*model.matrix(mod)[,15] + fixef(mod)[16]*model.matrix(mod)[,16] + fixef(mod)[17]*model.matrix(mod)[,17] + fixef(mod)[18]*model.matrix(mod)[,18] #+ fixef(mod)[19]*model.matrix(mod)[,19] + fixef(mod)[20]*model.matrix(mod)[,20] + fixef(mod)[21]*model.matrix(mod)[,21] + fixef(mod)[22]*model.matrix(mod)[,22] + fixef(mod)[23]*model.matrix(mod)[,23] + fixef(mod)[24]*model.matrix(mod)[,24] 

#Fixed <- fixef(mod)[2]*model.matrix(mod)[,2] 

VarF <- var(Fixed)

# R2GLMM(m) - marginal R2GLMM
# Equ. 26, 29 and 30
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
VarF/(VarF + VarCorr(mod)$Mno[1] + VarCorr(mod)$Study[1] + attr(VarCorr(mod), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. XXX, XXX
# my note: for variable slopes and intercepts, we want two parameters per random effect (intercept and logS). I don't think we indclude a correlation term.
# the output of VarCorr(mod)$Mno is variances/covariances, then stdev. Some useful code from when I didn't understand this: 
# print(summary(mod)$varcor) -> res
# as.numeric(res[,3]) -> res1
# attr(VarCorr(mod)$Mno, "stddev") would also work; except that gives STdevs, and we want variances, so stdevs squared.
# (VarF + res1[1]^2 + res1[2]^2 + res1[3]^2 + res1[4]^2)/(VarF + res1[1]^2 + res1[2]^2 + res1[3]^2 + res1[4]^2 + res1[5]^2)

(VarF + VarCorr(mod)$Mno[1] + VarCorr(mod)$Study[1])/(VarF + VarCorr(mod)$Mno[1] + VarCorr(mod)$Study[1] + attr(VarCorr(mod), "sc")^2)


## why R2GLMM(m) so much lower than what Brad gets for averaging individual fits?







##1. Preliminary examination of intercepts and slopes for each reported SST response
## following Pinheiro and Bates Ch 4
SST1.1<-as.data.frame(cbind(SST2[,1], SST2[,54], SST2[,51]))
#SST1.2<-SST1.1[-which(SST1.1$Mno==826),]
names(SST1.1)<-c('Mno', 'logY', 'logS')
SST1<-groupedData(logY~logS|Mno, SST1.1, order.groups=TRUE)
test<-lmList(logY~logS|Mno, data=SST1, na.action=na.omit)

pairs(test, id=0.01) 
head(coef(test))
mean(coef(test)[,2], na.rm=TRUE)

length(coef(test)[,2]) #230

plot(coef(test)[,1], coef(test)[,2], xlab='intercept', ylab='slope')
hist(coef(test)[,2], breaks=40, xlab='slope', main = 'Slopes from lmList individual model fitting')
abline(v=mean(coef(test)[,2]), col=2)

plot(intervals(test)) #about half the Mnos here have 2 richness levels, and i think that's preventing intervals from working. so i will make another file (above, SST3) with only studies that have >2 richness levels and do the whole anlyais on that, then rerun with teh full dataset.

### high overlap in confidence intervals suggests no need for random effects. We have large need for random intercepts, and possible need for random slopes.  





#estimating confidence intervals; this takes a while, and I"m not sure what the sigmas are.
mod<-modF1
confint(mod)
Computing profile confidence intervals ...
2.5 %    97.5 %
  .sig01       0.8620486 0.9904851
.sig02              NA        NA
.sig03              NA 0.0533406
.sig04       1.4965251 2.0741551
.sig05      -0.1867722 0.3124404
.sig06       0.2296943 0.3386019
.sigma       0.3227170 0.3496987
(Intercept)  4.5429167 5.3086167
logS         0.1661482 0.2986802

### ok, need to get the coefs for fixed and random effects for each observation (Mno | Study). the ranef extractor returns the conditional modes (which is our best estimate of the random effects)
str(ranef(modFM))
par(mfrow=c(1,2))
dotplot(ranef(modFM, condVar = TRUE))
qqmath(ranef(modFM, condVar = TRUE)) # figure this out

##### Producing experiment-level slope estimates #########
##########################################################

rand.cat <- data.frame(cbind(as.numeric(as.character(data$Mno)), as.numeric(as.character(data$Study))))
names(rand.cat) <- c('Mno', 'Study')
rand.cat1 <- ddply(rand.cat, .(Mno), summarize, mean(Study)) 

rand.cat <- ddply(data, .(Mno, Study, Sys1, TG, HigherT, restrt), summarize, mean(logY))
#rand.cat$Syst <- rand.cat$Sys1
names(rand.cat) <- c('Mno', 'Study', 'Syst','TG', 'HT', 'restrt', 'meanlogY')
Mno.coefs <- data.frame(coef(modBtrophic)$Mno)
Mno.coefs$Mno <- rownames(Mno.coefs)
#names(Mno.coefs)<-c('Int', 'LogS', 'System', 'HT', 'LogS.Sys1', 'LogS.HT', 'Mno1')
S <- cbind(rand.cat, Mno.coefs)

S$Sys.term <- ifelse(S$Syst == '1', S$logS.Sys1, 0)
S$TG.term<-ifelse(S$TG == '2', S$logS.TG12, 0)
S$TG.term<-ifelse(S$TG == '3', S$logS.TG13, S$TG.term)
S$TG.term<-ifelse(S$TG == '4', S$logS.TG14, S$TG.term)
#S$units <- ifelse(S$Units == 'density', S$logS.unit.types2density, 0)
#S$units <- ifelse(S$Units == 'perc.cover', S$logS.unit.types2perc.cover, S$units)
S$HT.term <- ifelse(S$HigherT == 'Y', S$logS.HT, 0)
#S$Res.term <- ifelse(S$restrt == 'incr', S$logS.restrtincr, 0)
#S$Res.term <- ifelse(S$restrt == 'red', S$logS.restrtred, S$Res.term)

#now add Study ranefs
St.ranefs <- data.frame(coef(modBtrophic)$Study)
St.ranefs$Study <- rownames(St.ranefs)
St.ranefs1 <- data.frame(St.ranefs$Study, St.ranefs$logS)
S2 <- merge(S, St.ranefs1, by.x = 'Study', by.y = 'St.ranefs.Study', all= FALSE)
S <- S2
b <- as.numeric(fixef(modBtrophic)[2])
S$slope <- S$logS + S$Sys.term + S$TG.term + S$HT.term + S$logS.log.Tscale. + (S$St.ranefs.logS - b)
hist(S$slope, breaks = 40, col = 'gray', freq = TRUE, main = 'Estimated b values', xlab = expression(beta[1][.ij] + beta[3][.ij] + mu[0][.i] + mu[0][.j] + mu[1][.i] + mu[1][.j]), xlim = c(-0.7, 1.4), ylim = c(0, 100), cex.lab = 1.2, axes = FALSE, ylab = 'Number of experiments') #
axis(1, at =c(-0.7, 0, 0.7, 1.4), lwd = 2, pos = 0)
axis(2, lwd = 2, pos = -0.7, las = 2)

summary(S$slope)  # ok, this is genreally working, though something is not quite right. the mean from this is 0.26, but the mean slope from the base model is 0.23. OK, the difference between these two is that the fixef for logS was counted twice - once in Mno and once in Study. that is now corrected by subtracting b in the final term above.
abline(v = mean(S$slope), lwd = 2)

par(new=T)
hist(S[(S$TG=='1'),]$slope, breaks = 20, xlim = c(-1, 1.5), ylim = c(0, 100), col = 3, axes = FALSE, xlab = '', ylab = '', main = '')

par(new=T)
hist(S[(S$TG=='4'),]$slope, breaks = 40, xlim = c(-1, 1.5), ylim = c(0, 100), col = 4, axes = FALSE, xlab = '', ylab = '', main = '')

par(new=T)
hist(S[(S$TG=='2'),]$slope, breaks = 20, xlim = c(-1, 1.5), ylim = c(0, 100), col = 5, axes = FALSE, xlab = '', ylab = '', main = '')

par(new=T)
hist(S[(S$TG=='3'),]$slope, breaks = 40, xlim = c(-1, 1.5), ylim = c(0, 100), col = 2, axes = FALSE, xlab = '', ylab = '', main = '')


#histograms with stacked bars
qplot(S$slope, binwidth = 0.025, fill = as.factor(S$TG))



summary(S$slope)

length(S[(S$TG=='3'),]$slope)
help(hist)


### calculating CIs for the estimates:
n<-length(S$slope)
est<- mean(S$slope)
se <- sd((S$slope/sqrt(n)))
in.95 <- est + qt(c(0.025, 0.975), n-1)*se
abline(v = in.95[1], lwd = 2, lty = 2)
abline(v = in.95[2], lwd = 2, lty = 2)



#######
#trophic group summary: what are the carnivore data?
#checking area x duration relationship
lm(log(SST2$MaxTscale)~log(SST2$Vol)) -> vol.mod
lm(log(SST2$MaxTscale)~log(SST2$Area)) -> area.mod
summary(vol.mod)
summary(area.mod)

par(mfrow = c(2,1))
plot(log(SST2$MaxTscale)~log(SST2$Vol))
abline(coef(vol.mod[1]),coef(vol.mod[2]))
text(2, 4, 'p < 0.001, m = -0.4')

plot(log(SST2$MaxTscale)~log(SST2$Area))
abline(coef(area.mod[1]),coef(area.mod[2]))
text(2, 4, 'p < 0.001, m = 0.2')

par(mfrow = c(2,1))
plot(log(SST2$MaxTscale)~log(SST2$SPscale))
lm(log(SST2$MaxTscale)~log(SST2$SPscale+1)) -> spcalemod
abline(coef(vol.mod[1]),coef(vol.mod[2]))
text(2, 4, 'p < 0.001, m = -0.4')

plot(log(SST2$MaxTscale)~log(SST2$Area))
abline(coef(area.mod[1]),coef(area.mod[2]))
text(2, 4, 'p < 0.001, m = 0.2')



#### plotting the power function
a <- 1
b <- .5
seq(0, 200, 0.1) -> x
pwer <- function(x) a*x^b
plot(pwer(x)~x, ylim = c(0, 120))