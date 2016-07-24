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

data <- read.csv("/Users/maryo/Documents/projects/BEF synthesis/BEF-scaling-MA/data/SST5.csv")
data <- SST5

library(lme4)
library(plyr)

mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

## create the Entry, Study and Expt columns for later matching with random effect values
rand.cat <- ddply(data, .(Entry, Study, ExptA, Sys1, TG1), summarize, mean(logY.rs))
names(rand.cat) <- c('Entry', 'Study', 'ExptA','System','TG1', 'meanlogY')
Entry.coefs <- data.frame(coef(mod4)$Entry) #this gives us the logS for entry. so we'd need to add only the expt and study level later. let's try that. AND, I can see that logSc here already includes the random effect
Entry.coefs$Entry <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by = 'Entry', all = FALSE)

b <- as.numeric(fixef(mod4)[2])
se <- function(x) sd(x)/sqrt(length(x))

## constructing predicted b values (slopes) from fixed effects for each group.
S$Trt <- paste(S$System, S$TG1)
S$Trt.term <- ifelse(S$Trt == 'A 1', S$logSc, 0)
S$Trt.term <- ifelse(S$Trt == 'A 2', (S$logSc + S$logSc.TG12), S$Trt.term)
S$Trt.term <- ifelse(S$Trt == 'A 4', (S$logSc + S$logSc.TG14), S$Trt.term)
S$Trt.term <- ifelse(S$Trt == 'T 1', (S$logSc + S$logSc.Sys1T), S$Trt.term)
S$Trt.term <- ifelse(S$Trt == 'T 4', (S$logSc + S$logSc.TG14 + S$logSc.Sys1T.TG14), S$Trt.term)

## producing Treatment mean coefficients
b.sum <- ddply(S, .(Trt), summarize, mean(Trt.term))
b.sum2 <- ddply(S, .(Trt), summarize, se(Trt.term))
b.sums <- merge(b.sum, b.sum2, by = 'Trt')
names(b.sums) <- c('group', 'est', 'se')

b.sums$group <- c("Aq. Primary Prod.", "Aq. Herbivores", "Aq. Detritovores", "Terr. Prim. Prod", "Terr. Detrit.")

### two-paneled figure

pdf(file = "figure 1B.pdf", width = 4, height = 4)
par(
  family = "serif",  
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  #fin = c(7,5), pty = "m",
  mar = c(5,10,4,0),  # Inner margins are set through a little trial and error.
  mfcol = c(1,1)
)

#TOP PANEL: SLOPES
par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(-0.2, 0.6),                        	
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
  text(0.6, length(b.sums[,1])+ .2, 'B', cex = 1.2)
}

# add axes and labels
axis(side = 1)                                                                                         
#abline(v = 0, lty = 3, col = "grey40")                                     
               
#abline(h = 3.5, lty = 3, col = 'grey40')     
#abline(h = 5.5, lty = 3, col = 'grey40')  

mtext(side = 1, "Slope coefficients", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          

dev.off()

#### repeat for intercepts
## constructing predicted b values (slopes) from fixed effects for each group.
#S$Trt <- paste(S$Syst, S$TG1)
S$ITrt.term <- ifelse(S$Trt == 'A 1', S$X.Intercept., 0)
S$ITrt.term <- ifelse(S$Trt == 'A 2', (S$X.Intercept. + S$TG12), S$ITrt.term)
S$ITrt.term <- ifelse(S$Trt == 'A 4', (S$X.Intercept. + S$TG14), S$ITrt.term)
S$ITrt.term <- ifelse(S$Trt == 'T 1', (S$X.Intercept. + S$Sys1T), S$ITrt.term)
S$ITrt.term <- ifelse(S$Trt == 'T 4', (S$X.Intercept. + S$TG14 + S$Sys1T.TG14), S$ITrt.term)

## producing Treatment mean coefficients
I.sum <- ddply(S, .(Trt), summarize, mean(ITrt.term))
I.sum2 <- ddply(S, .(Trt), summarize, se(ITrt.term))
I.sums <- merge(I.sum, I.sum2, by = 'Trt')
names(I.sums) <- c('group', 'est', 'se')

I.sums$group <- c("Aq. Primary Prod.", "Aq. Herbivores", "Aq. Detritovores", "Terr. Prim. Prod", "Terr. Detrit.")



### INTERCEPTS
pdf(file = "figure 1B ints.pdf", width = 4, height = 4)

par(
  family = "serif",  
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  #fin = c(7,5), pty = "m",
  mar = c(5,10,4,0),  # Inner margins are set through a little trial and error.
  mfcol = c(1,1)
)

#TOP PANEL: SLOPES
par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(2, 6),                        	
     ylim = c(.7, length(I.sums[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests.I <- (I.sums[,2])
ses.I <- (I.sums[,3])
var.names <- (I.sums[,1])

a <- 0
for (i in 1:length(I.sums[,1])) {                                            
  points(ests.I[i], i+a, pch = 19, cex = 1.2, col = 1) 
  lines(c(ests.I[i] + 1.96*ses.I[i], ests.I[i] - 1.96*ses.I[i]), c(i+a, i+a), col = 'gray80', lwd = 3)
  text(1.5, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names
  text(0.6, length(I.sums[,1])+ .2, '', cex = 1.2)
}

# add axes and labels
axis(side = 1)                                                                                         
#abline(v = 0, lty = 3, col = "grey40")                                     
abline(v = 0, lty = 1, col = 'grey60')               
#abline(h = 3.5, lty = 3, col = 'grey40')     
#abline(h = 5.5, lty = 3, col = 'grey40')  

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



