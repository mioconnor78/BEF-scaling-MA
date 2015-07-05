#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure 2
### Jan 06 2015; Author: Mary O'Connor
### updated Jan 26 for modBasic and modBtrophic
#######################################################################

## a four panel figure: 2 cols for w/ and w/o preds, and 2 rows for intercept and slope coefs
##############################################################################################

## MO: I know how painful and inefficient this code is! I would love it (and learn from it) if someone felt like presenting an alternative approach to generating this figure. I know there are better ways.

## code for averaged model with 3-way interaction here. Below (waaaay down) is code for previous best model that didn't have this interaction.


### CODE FOR THE AVERAGE MODEL (4 MODELS)

#layout(matrix(c(1,2,3,4)), 1, 1, byrow = FALSE)

# Figure 2A: SST5 slopes
# create file for model averaged estimates
## July 5: I'm going to skip the model averaging for now, b/c no delta AIC < 2.
#estimates <- as.data.frame(m.avg[3])
#rownames(estimates) <- c('Intercept', 'ln(S)', 'ln(Tg)', 'Ecosystem', 'Herbivore', 'Detritivore', 'ln(S)*ln(Tg)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore','Ecosyst*Detritivore', 'ln(S)*Syst*Detrit', '+Resource', '-Resource', 'ln(maxTg)', 'ln(S) * +Resource', 'ln(S)* -Resource', 'ln(S)* ln(maxTg)', 'Density','Percent Cover','LabField','ln(Smax)','ln(S)*Density','ln(S)*Percent Cover','ln(S)*LabField', 'ln(S)*ln(Smax)' )
#colnames(estimates) <- c('est', 'se', 'adjse', 'lCI', 'uCI')
#estimates$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'I', 'S', 'I', 'I','I', 'S', 'S', 'S','I', 'I', 'I','I', 'S', 'S', 'S','S')
#est.sl <- estimates[estimates$slint == 'S',]
#est.int <- estimates[estimates$slint == 'I',]

# create data for best model estimates
est.B <- as.data.frame(as.numeric(round(fixef(mod4.1),3)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(mod4.1))),3))
names(est.B) <- c('est', 'se')
rownames(est.B) <- c('Intercept', 'ln(Tg)', 'ln(S)', 'Ecosystem', 'Herbivore', 'Detritivore', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore','Ecosyst*Detrit', 'ln(S)*Syst*Detrit')
est.B$slint <- c('I', 'I', 'S',  'I', 'I', 'I', 'S', 'S', 'S', 'I', 'S')
est.B.sl <- est.B[est.B$slint == 'S',]
est.B.int <- est.B[est.B$slint == 'I',]


# create data for basic model estimates  #needs review
est.Ba <- as.data.frame(as.numeric(round(fixef(mod3),3)))
est.Ba$se <- as.numeric(round(sqrt(diag(vcov(mod3))),3))
names(est.Ba) <- c('est', 'se')
est.Ba <- rbind(est.Ba[1,], c('',''), est.Ba[2,], c('',''), c('',''), c('',''), c('',''),  c('',''), c('',''), c('',''), c('','')) 
rownames(est.Ba) <- c('Intercept', 'ln(Tg)', 'ln(S)', 'Ecosystem', 'Herbivore', 'Detritivore', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore','Ecosyst*Detrit', 'ln(S)*Syst*Detrit')
est.Ba$slint <- c('I', 'I', 'S',  'I', 'I', 'I', 'S', 'S', 'S', 'I', 'S')
est.Ba.sl <- est.Ba[est.Ba$slint == 'S',]
est.Ba.int <- est.Ba[est.Ba$slint == 'I',]

### two-paneled figure

pdf(file = "figure 1B.pdf", width = 7.5, height = 4)
par(
  family = "serif",  
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  #fin = c(7,5), pty = "m",
  mar = c(5,10,4,0),  # Inner margins are set through a little trial and error.
  mfcol = c(1,2)
)

#TOP PANEL: SLOPES
par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(-1.0, 0.6),                        	
     ylim = c(.7, length(est.B.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
#est <- as.numeric(est.sl[,1]) 
#se <- as.numeric(est.sl[,2] )                                         
ests.B <- as.numeric(est.B.sl[,1])
ses.B <- as.numeric(est.B.sl[,2])
ests.Ba <- as.numeric(est.Ba.sl[,1])
ses.Ba <- as.numeric(est.Ba.sl[,2])
var.names<-rownames(est.B.sl)
var.namesi<-rownames(est.B.int)

b <- 0
for (i in 1:length(ests.B)) {                                            
  #points(est[i], i, pch = 19, cex = 1.2)                              
  #lines(c(est[i] + 1.96*se[i], est[i] - 1.96*se[i]), c(i, i), lwd = 2)         # add 95% CIs
  points(ests.B[i], i+b, pch = 19, cex = 1.2, col = 1) 
  lines(c(ests.B[i] + 1.96*ses.B[i], ests.B[i] - 1.96*ses.B[i]), c(i+b, i+b), col = 1, lwd = 2)
  points(ests.Ba[i], i+2*b, pch = 19, cex = 1.2, col = 'gray50') 
  lines(c(ests.Ba[i] + 1.96*ses.Ba[i], ests.Ba[i] - 1.96*ses.Ba[i]), c(i+2*b, i+2*b), col = 'gray50', lwd = 2)
  text(-1.2, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names
  text(0.45, length(est.B.sl[,1])+ .2, 'B', cex = 1.2)
}

# add axes and labels
axis(side = 1)                                                                                         
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Slope coefficients", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          

### INTERCEPTS
par(mar=c(5,8,4,4))  #pin = c(2.3, 3.5)), but this doesn't seem to work with mar
plot(NULL,                                
     xlim = c(-6, 6),                          
     ylim = c(.7, length(est.B.int[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
#est <- as.numeric(est.int[,1]) 
#se <- as.numeric(est.int[,2] )                                         
ests.B <- as.numeric(est.B.int[,1])
ses.B <- as.numeric(est.B.int[,2])
ests.Ba <- as.numeric(est.Ba.int[,1])
ses.Ba <- as.numeric(est.Ba.int[,2])
var.names<-rownames(est.B.int)

b <- 0
for (i in 1:length(ests.B)) {                                            
  #points(est[i], i, pch = 19, cex = 1.2)                              
  #lines(c(est[i] + 1.96*se[i], est[i] - 1.96*se[i]), c(i, i), lwd = 2)  
  points(ests.B[i], i+b, pch = 19, cex = 1.2, col = 1) 
  lines(c(ests.B[i] + 1.96*ses.B[i], ests.B[i] - 1.96*ses.B[i]), c(i+b, i+b), col = 1, lwd = 2)
  lines(c(ests.Ba[i] + 1.96*ses.Ba[i], ests.Ba[i] - 1.96*ses.Ba[i]), c(i+2*b, i+2*b), col = 'gray50', lwd = 2)
  points(ests.Ba[i], i+2*b, pch = 19, cex = 1.2, col = 'gray50')   # add 95% CIs
  text(-7, i, adj = c(1,0), var.namesi[i], xpd = T, cex = .8)        # add the variable names
  text(5.5, length(est.B.int[,1])+ 0.2, 'C', cex = 1.2)
}

# add axes and labels
axis(side = 1, at = c(-6, 0, 6))
#axis(side = 2, pos = -2)
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Intercept coefficients", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                    

dev.off()












###  OLD CODE FOR AN AVERAGED MODEL WITHOUT THE 3-WAY INTERACTION ####

#layout(matrix(c(1,2,3,4)), 1, 1, byrow = FALSE)

# Figure 2A: SST4 slopes
# create file for model averaged estimates
estimates <- as.data.frame(m.avg[2])
rownames(estimates) <- c('Intercept', 'ln(S)', 'ln(Tg)', 'Ecosystem', 'Herbivore', 'Detritivore', 'ln(S)*ln(Tg)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore', '+Resource', '-Resource','ln(maxTg)', 'ln(S) * +Resource', 'ln(S)* -Resource', 'ln(S)* ln(maxTg)')
colnames(estimates) <- c('est', 'se', 'adjse', 'lCI', 'uCI')
estimates$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'I', 'I','I', 'S', 'S', 'S')
est.sl <- estimates[estimates$slint == 'S',]
est.int <- estimates[estimates$slint == 'I',]

# create data for best model estimates
est.B <- as.data.frame(as.numeric(round(fixef(modBtrophic),3)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(modBtrophic))),3))
names(est.B) <- c('est', 'se')
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
rownames(est.B) <- c('Intercept', 'ln(S)', 'ln(Tg)', 'Ecosystem', 'Herbivore', 'Detritivore', 'ln(S)*ln(Tg)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore', '+Resource', '-Resource','ln(maxTg)', 'ln(S) * +Resource', 'ln(S)* -Resource', 'ln(S)* ln(maxTg)')
est.B$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S','I', 'I','I','S', 'S', 'S')
est.B.sl <- est.B[est.B$slint == 'S',]
est.B.int <- est.B[est.B$slint == 'I',]


# create data for basic model estimates
est.Ba <- as.data.frame(as.numeric(round(fixef(modBasic),3)))
est.Ba$se <- as.numeric(round(sqrt(diag(vcov(modBasic))),3))
names(est.Ba) <- c('est', 'se')
est.Ba <- rbind(est.Ba[1:3,], c('',''), c('',''), c('',''), est.Ba[4:nrow(est.Ba),], c('',''), c('',''), c('',''), c('',''), c('',''), c('',''), c('',''), c('',''), c('','')) 
rownames(est.Ba) <- c('Intercept', 'ln(S)', 'ln(Tg)', 'Ecosystem', 'Herbivore', 'Detritivore', 'ln(S)*ln(Tg)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Detritivore',  '+Resource', '-Resource', 'ln(maxTg)', 'ln(S) * +Resource', 'ln(S)* -Resource', 'ln(S)* ln(maxTg)')
est.Ba$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S','I', 'I','I','S', 'S', 'S')
est.Ba.sl <- est.Ba[est.Ba$slint == 'S',]
est.Ba.int <- est.Ba[est.Ba$slint == 'I',]

### two-paneled figure

par(
  family = "serif",  
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  #mar = c(5,8,4,2),  # Inner margins are set through a little trial and error.
  mfcol = c(1,2)
)

#TOP PANEL: SLOPES

par(mar=(c(5,9,4,0)))
plot(NULL,                                
     xlim = c(-0.3, 0.4),                        	
     ylim = c(.7, length(est.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
est <- as.numeric(est.sl[,1]) 
se <- as.numeric(est.sl[,2] )                                         
ests.B <- as.numeric(est.B.sl[,1])
ses.B <- as.numeric(est.B.sl[,2])
ests.Ba <- as.numeric(est.Ba.sl[,1])
ses.Ba <- as.numeric(est.Ba.sl[,2])
var.names<-rownames(est.B.sl)
var.namesi<-rownames(est.B.int)

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = 1.2)                              
  points(ests.B[i], i+b, pch = 19, cex = 1.2, col = 'gray50') 
  points(ests.Ba[i], i+2*b, pch = 19, cex = 1.2, col = 'gray75') 
  lines(c(est[i] + 1.96*se[i], est[i] - 1.96*se[i]), c(i, i), lwd = 2)         # add 95% CIs
  lines(c(ests.B[i] + 1.96*ses.B[i], ests.B[i] - 1.96*ses.B[i]), c(i+b, i+b), col = 'gray50', lwd = 2)
  lines(c(ests.Ba[i] + 1.96*ses.Ba[i], ests.Ba[i] - 1.96*ses.Ba[i]), c(i+2*b, i+2*b), col = 'gray75', lwd = 2)
  text(-0.4, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names
  text(0.35, length(est.B.sl[,1]), 'A', cex = 1.5)
}

# add axes and labels
axis(side = 1)                                                                                         
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Slope coefficients", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          


par(mar=(c(5,5,4,4)))
plot(NULL,                                
     xlim = c(-3, 6),                          
     ylim = c(.7, length(est.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
est <- as.numeric(est.int[,1]) 
se <- as.numeric(est.int[,2] )                                         
ests.B <- as.numeric(est.B.int[,1])
ses.B <- as.numeric(est.B.int[,2])
ests.Ba <- as.numeric(est.Ba.int[,1])
ses.Ba <- as.numeric(est.Ba.int[,2])
var.names<-rownames(est.B.int)

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = 1.2)                              
  points(ests.B[i], i+b, pch = 19, cex = 1.2, col = 'gray50') 
  points(ests.Ba[i], i+2*b, pch = 19, cex = 1.2, col = 'gray75') 
  lines(c(est[i] + 1.96*se[i], est[i] - 1.96*se[i]), c(i, i), lwd = 2)         # add 95% CIs
  lines(c(ests.B[i] + 1.96*ses.B[i], ests.B[i] - 1.96*ses.B[i]), c(i+b, i+b), col = 'gray50', lwd = 2)
  lines(c(ests.Ba[i] + 1.96*ses.Ba[i], ests.Ba[i] - 1.96*ses.Ba[i]), c(i+2*b, i+2*b), col = 'gray75', lwd = 2)
  text(-4, i, adj = c(1,0), var.namesi[i], xpd = T, cex = .8)        # add the variable names
  text(5.5, length(est.B.int[,1]), 'B', cex = 1.5)
}

# add axes and labels
axis(side = 1, at = c(-2, 0, 2, 4, 6))
#axis(side = 2, pos = -2)
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Intercept coefficients", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                    



