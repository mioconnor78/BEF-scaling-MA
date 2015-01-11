#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure 2
### Jan 06 2015; Author: Mary O'Connor
#######################################################################

## a four panel figure: 2 cols for w/ and w/o preds, and 2 rows for intercept and slope coefs
##############################################################################################

## MO: I know how painful and inefficient this code is! I would love it (and learn from it) if someone felt like presenting an alternative approach to generating this figure. I know there are better ways.


par(
  family = "serif",  
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  #mar = c(5,8,4,2),  # Inner margins are set through a little trial and error.
  mfcol = c(2,2)
)

#layout(matrix(c(1,2,3,4)), 1, 1, byrow = FALSE)

# Figure 2A: SST4 slopes
estimates <- as.data.frame(m.avg.4[3])
estimates$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'S', 'S', 'I', 'I','S', 'S')
rownames(estimates) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Predator', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator','ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
est.sl <- estimates[estimates$slint == 'S',]
est.int <- estimates[estimates$slint == 'I',]

est.B <- as.data.frame(as.numeric(round(fixef(modBtrophic4),3)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(modBtrophic4))),3))
names(est.B) <- c('est', 'se')
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
row.names(est.B) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Predator', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator','ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
est.B$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'S', 'S', 'I', 'I','S', 'S')

est.B.sl <- est.B[est.B$slint == 'S',]
est.B.int <- est.B[est.B$slint == 'I',]

#par(mar=(c(5,7.5,4,0)))
par(mar=(c(5,9.5,4,0)))
plot(NULL,                                
     xlim = c(-0.4, 0.6),                      		
     ylim = c(.7, length(est.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
est <- as.numeric(est.sl[,1]) 
se <- as.numeric(est.sl[,2] )                                         
ests.B <- as.numeric(est.B.sl[,1])
ses.B <- as.numeric(est.B.sl[,2])
var.names<-rownames(est.B.sl)

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = .75)                              
  points(ests.B[i], i+b, pch = 19, cex = .75, col = 'gray60') 
  lines(c(est[i] + 1.94*se[i], est[i] - 1.94*se[i]), c(i, i))         # add 90% CIs
  lines(c(ests.B[i] + 1.94*ses.B[i], ests.B[i] - 1.94*ses.B[i]), c(i+b, i+b), col = 'gray60')         
  text(-.5, i, adj = c(1,0), var.names[i], xpd = T, cex = 1)        # add the variable names
  text(0.55, length(est.B.sl[,1]), 'A', cex = 1.5)
}

# add axes and labels
axis(side = 1)                                                                                         
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Slope coefficients", line = 3)                                              
mtext(side = 3, "Predator Studies Included", line = 1)   # add title
box()                                          

## Panel C: SST4 Ints ###
par(mar=(c(8,9.5,1,0)))
plot(NULL,                                
     xlim = c(-3, 4),                        	
     ylim = c(.7, length(est.int[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
est <- as.numeric(est.int[,1]) 
se <- as.numeric(est.int[,2] )                                         
ests.B <- as.numeric(est.B.int[,1])
ses.B <- as.numeric(est.B.int[,2])
var.names<-rownames(est.B.int)

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = .75)                              
  points(ests.B[i], i+b, pch = 19, cex = .75, col = 'gray60') 
  lines(c(est[i] + 1.94*se[i], est[i] - 1.94*se[i]), c(i, i))         # add 90% CIs
  lines(c(ests.B[i] + 1.94*ses.B[i], ests.B[i] - 1.94*ses.B[i]), c(i+b, i+b), col = 'gray60')         
  text(-3.7, i, adj = c(1,0), var.names[i], xpd = T, cex = 1)        # add the variable names
  text(3.5, length(est.B.int[,1]), 'C', cex = 1.5)
}

# add axes and labels
axis(side = 1)                                                                                      
abline(v = 0, lty = 3, col = "grey40")                                                                    
mtext(side = 1, "Intercept Coefficients", line = 3)                                              
#mtext(side = 3, "C. With Predators", line = 1)   # add title
box()                    


## Panel B
# SST5 slopes
estimates <- as.data.frame(m.avg.5[3])
estimates <- rbind(estimates[1:5,], c('','','','','',''), estimates[6:nrow(estimates),])  #these rows only needed for SST5 (no carnivores)
estimates <- rbind(estimates[1:11,], c('','','','','',''), estimates[12:nrow(estimates),])
estimates$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'S', 'S', 'I', 'I','S', 'S')
rownames(estimates) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Predator', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator','ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
est.sl <- estimates[estimates$slint == 'S',]
est.int <- estimates[estimates$slint == 'I',]

est.B <- as.data.frame(as.numeric(round(fixef(modBtrophic5),3)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(modBtrophic5))),3))
names(est.B) <- c('est', 'se')

#the parameters in estimates and est.B have to line up
rownames(est.B)
rownames(estimates)

est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B, c('',''))
est.B <- rbind(est.B[1:5,], c('','','','','',''), est.B[6:nrow(est.B),]) #these rows only needed for SST5 (no carnivores)
est.B <- rbind(est.B[1:11,], c('','','','','',''), est.B[12:nrow(est.B),])
row.names(est.B) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Predator', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator','ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
est.B$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'S', 'S', 'I', 'I','S', 'S')
est.B.sl <- est.B[est.B$slint == 'S',]
est.B.int <- est.B[est.B$slint == 'I',]


par(mar=(c(5,5.5,4,4)))
plot(NULL,                                
     xlim = c(-0.4, 0.6),                        	
     ylim = c(.7, length(est.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
est <- as.numeric(est.sl[,1]) 
se <- as.numeric(est.sl[,2] )                                         
ests.B <- as.numeric(est.B.sl[,1])
ses.B <- as.numeric(est.B.sl[,2])
var.names<-rownames(est.B.sl)

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = .75)                              
  points(ests.B[i], i+b, pch = 19, cex = .75, col = 'gray60') 
  lines(c(est[i] + 1.94*se[i], est[i] - 1.94*se[i]), c(i, i))         # add 90% CIs
  lines(c(ests.B[i] + 1.94*ses.B[i], ests.B[i] - 1.94*ses.B[i]), c(i+b, i+b), col = 'gray60')         
  #text(-.8, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names
  text(0.55, length(est.sl[,1]), 'B', cex = 1.5)
}

# add axes and labels
axis(side = 1)                                                                                      
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Slope coefficients", line = 3)                                             
mtext(side = 3, "Predator Studies Excluded", line = 1)   # add title
box()                            


### Panel D: SST5 intercepts
par(mar=(c(8,5.5,1,4)))
plot(NULL,                                
     xlim = c(-3, 4),                          
     ylim = c(.7, length(est.B.int[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
est <- as.numeric(est.int[,1]) 
se <- as.numeric(est.int[,2] )                                         
ests.B <- as.numeric(est.B.int[,1])
ses.B <- as.numeric(est.B.int[,2])
var.names<-rownames(est.B.int)

b <- 0.2
for (i in 1:length(est)) {                                            
  points(est[i], i, pch = 19, cex = .75)                              
  points(ests.B[i], i+b, pch = 19, cex = .75, col = 'gray60') 
  lines(c(est[i] + 1.94*se[i], est[i] - 1.94*se[i]), c(i, i))         # add 90% CIs
  lines(c(ests.B[i] + 1.94*ses.B[i], ests.B[i] - 1.94*ses.B[i]), c(i+b, i+b), col = 'gray60')         
  #text(-3.6, i, adj = c(1,0), var.names[i], xpd = T, cex = 1)     
  text(3.5, length(est.B.int[,1]), 'D', cex = 1.5)
}

# add axes and labels
axis(side = 1)                                                                                        
abline(v = 0, lty = 3, col = "grey40")                                                                 
mtext(side = 1, "Intercept coefficients", line = 3)                                             
#mtext(side = 3, "D. Without Predators", line = 1)   
box()                    


###### END OF PLOT #####






## general data preparation, though it is embedded above so this code is extra
## plotting from http://www.carlislerainey.com/Blog_Files/Blog_CoefficientPlots.R

estimates <- as.data.frame(m.avg[3])
estimates <- rbind(estimates[1:5,], c('','','','','',''), estimates[6:nrow(estimates),])  #these rows only needed for SST5 (no carnivores)
estimates <- rbind(estimates[1:11,], c('','','','','',''), estimates[12:nrow(estimates),])
estimates$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'S', 'S', 'I', 'I','S', 'S')
rownames(estimates) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Predator', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator','ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
est.sl <- estimates[estimates$slint == 'S',]
est.int <- estimates[estimates$slint == 'I',]


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
est.B <- rbind(est.B[1:5,], c('','','','','',''), est.B[6:nrow(est.B),]) #these rows only needed for SST5 (no carnivores)
est.B <- rbind(est.B[1:11,], c('','','','','',''), est.B[12:nrow(est.B),])
row.names(est.B) <- c('Intercept', 'ln(S)', 'ln(Duration)', 'Ecosystem', 'Herbivore', 'Predator', 'Detritivore', '+Consumer', 'ln(S)*ln(Duration)', 'ln(S)*Ecosystem', 'ln(S)*Herbivore', 'ln(S)*Predator','ln(S)*Detritivore', 'ln(S)* +Consumer', '+Resource', '-Resource', 'ln(S) * +Resource', 'ln(S)* -Resource')
est.B$slint <- c('I', 'S', 'I', 'I', 'I', 'I', 'I', 'I', 'S', 'S', 'S', 'S', 'S', 'S', 'I', 'I','S', 'S')
#est.B <- as.data.frame(est.B)

est.B.sl <- est.B[est.B$slint == 'S',]
est.B.int <- est.B[est.B$slint == 'I',]


