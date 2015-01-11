#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure 2
### Jan 06 2015; Author: Mary O'Connor
#######################################################################

## requires model results as input
## plotting from http://www.carlislerainey.com/Blog_Files/Blog_CoefficientPlots.R
## ok, we want a 2-panel figure, one panel for the w/ the full best model and one for model avg.coefs. one set can be in grey. let's try that. 

# a full figure with intercept coefs
estimates <- as.data.frame(m.avg[3])
var.names<-rownames(estimates)
var.names<-c('Intercept', 'logS', 'Ecosystem - Aq', 'TG: Herbivore', 'TG: Carnivore', 'TG: Detritovore', 'Consumer present', '+ Resource', '- Resource', 'log(Time in gen)', 'logS*Ecosystem', 'logS*Herbivore', 'logS*Carnivore', 'logS*Detritovore', 'logS*Consumer pres', 'logS * + Resource', 'logS* - Resource', 'logS*log(Time in gen)', 'units Density', 'units % Cover', 'log(Smax)', 'logS*Density', 'logS*mass', 'logS*log(Smax)')

est.B <- as.data.frame(as.numeric(round(fixef(modBtrophic),3)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(modBtrophic))),3))
names(est.B) <- c('est', 'se')

#the parameters in estimates and est.B have to line up
rownames(est.B)
rownames(estimates)

est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
row.names(est.B) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
#est.B <- as.data.frame(est.B)
var.names<-rownames(est.B)

# set the graphical parameters
par(
  family = "serif",  # I don't plot in anything but serif
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  mar = c(5,10,4,2)  # Inner margins are set through a little trial and error.
)

# create an empty plot for total customization
plot(NULL,                                  # create empty plot
     xlim = c(-3, 4),                      		# set xlim by guessing
     ylim = c(.7, length(estimates[,1]) + .3), 	# set ylim by the number of variables
     axes = F, xlab = NA, ylab = NA)       		# turn off axes and labels

# add the data
est <- estimates[,1] #[-1] # conveniently store the estimates (minus the constant)
se <- estimates[,2]                                        		# conveniently store the std. errors (minus the constant)
ests.B <- as.numeric(est.B[,1])
ses.B <- as.numeric(est.B[,2])

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = .75)                              
  points(ests.B[i], i+b, pch = 19, cex = .75, col = 'gray60') 
  lines(c(est[i] + 1.94*se[i], est[i] - 1.94*se[i]), c(i, i))         # add 90% CIs
  lines(c(ests.B[i] + 1.94*ses.B[i], ests.B[i] - 1.94*ses.B[i]), c(i+b, i+b), col = 'gray60')         
  #lines(c(estimates$avg.model.Upper.CI[i], estimates$avg.model.Lower.CI[i]), c(i, i))
  text(-3.6, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names

}

# add axes and labels
axis(side = 1)                                                                                          # add bottom axis
abline(v = 0, lty = 3, col = "grey40")                                                                    # add verticle line
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


