### estimating proportion of function remaining after species loss

## proportion remaining:

#Y.R <- 1 - (S.2^b)/(S.1^b)

# do some algebra, and: 

b <- 0.26 # scaling exponent
S.1 <- seq(1, 100, 1) #species richness before loss
#S.2 number of species lost
#Y.R proportion biomass remaining after species loss

S.2.func = function(Y, S.1) exp((1/b) * ( log(Y) + b*log(S.1) ))

### so what if we plot a histogram of expected biomass change due to change in diversity. This would be acheived by rearranging the S.2 function to solve for Y2, for change in richness and the starting richness. 
Y.func <- function(S.1, S.2) (S.2/S.1)^b  # Y.postloss / Y.preloss = proportion of biomass remaining


#how many species do you have to lose to lose 10% of function?
pdf(file = "figure RB.pdf", width = 6, height = 6)
plot(S.2.func(1, S.1) ~ S.1, xlim = c(1,100), ylim = c(0,100), pch = '', xlab = 'starting species richness', ylab = 'Species required for min biomass (YR) after sp loss', cex.lab = 1.5)
lines(S.2.func(1, S.1) ~ S.1, col=1, lwd = 2)
lines(S.2.func(0.9, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(0.8, S.1) ~ S.1, col='gray40', lwd = 2)
lines(S.2.func(0.5, S.1) ~ S.1, col='gray60', lwd = 2)
lines(S.2.func(0.1, S.1) ~ S.1, col='gray70', lwd = 2)
lines(S.2.func(0, S.1) ~ S.1, col='gray80', lwd = 2)
legend(-5, 100, c('YR = 1', 'YR = 0.9', 'YR = 0.8', 'YR = 0.5', 'YR = 0.1', 'YR = 0'), pch=19, col = c('black', 'gray20', 'gray40',  'gray60', 'gray70','gray90'), bty = 'n')
abline(v = 100, lty = 2, lwd = 1)
abline(v = 20, lty = 2, lwd = 1)
dev.off()

#ok, so with Y.R set at 0, we see a 1:1 line. meaning, to maintain 100% function you can lose no species. 
# Y.R = 0.05, you can have 80 species and acheive the 95% of the function of 100 specie.
# so, if Y.R = 0.5, we see that 50% function is acheived by something like 10 species; in other words, if you start with 100 species, you can lose 90 and maintain 50% of function. If you start with 40 species, you can end up with as few as 4 (or so) and still maintain 50% function. 
# if Y.R = 0.9, i.e., the post-loss community retains 10% of function, then you can lose 


### bring in vellend data

data=read.csv("./vellend.csv")
data$dSR <- data$SR_Year1_CT - data$SR_Year2_CT
hist(data$dSR, breaks = 40)
hist(data$SR_Year1_CT, breaks = 40)
plot(data$dSR ~ data$SR_Year1_CT)

pdf(file = "figure RB.pdf", width = 6, height = 6)
plot(data$SR_Year2_CT ~ data$SR_Year1_CT, ylim = c(0,100), xlim = c(0,100), xlab = 'Species Richness before change', ylab = 'Species Richness after change')
abline(a = 1, b = 1, lwd = 2)
lines(S.2.func(0.9, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(0.8, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(0.5, S.1) ~ S.1, col='red', lwd = 2)
lines(S.2.func(1.1, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(1.2, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(1.5, S.1) ~ S.1, col='red', lwd = 2)
legend(-5, 105, c('YR = 1', 'YR = 0.9', 'YR = 0.8', 'YR = 0.5', 'YR = 1.1', 'YR = 1.2', 'YR = 1.5'), pch=19, col = c('black', 'gray20', 'gray40',  'gray60', 'gray20', 'gray40',  'gray60'), bty = 'n')
dev.off()




pdf(file = "figure deltaY Vellend.pdf", width = 4, height = 4)
hist(Y.func(data$SR_Year2_CT, data$SR_Year1_CT), col = 'gray40', breaks = 40, main = 'Change in function predicted for Vellend plots', xlab = 'Expected Biomass Change (Y.2/Y.1)', las = 1)
abline(v = 0.9, lwd = 2)
abline(v = 1.1, lwd = 2)
dev.off()


# this shows that it's not the mangitude of change, but the combination of magnitude and where your'e starting from.
plot(Y.func(data$SR_Year2_CT, data$SR_Year1_CT) ~ I(data$SR_Year2_CT - data$SR_Year1_CT), ylim = c(0,2), xlim = c(-20,50), xlab = 'Species Richness change', ylab = 'proportion biomass change resulting from richness change')

plot(Y.func(data$SR_Year2_CT, data$SR_Year1_CT) ~ data$SR_Year1_CT, ylim = c(0,2), xlim = c(0,100), xlab = 'Species Richness before change', ylab = 'proportion biomass change resulting from richness change')

abline(a = 1, b = 1, lwd = 2)
lines(S.2.func(0.9, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(0.8, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(0.5, S.1) ~ S.1, col='red', lwd = 2)
lines(S.2.func(1.1, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(1.2, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(1.5, S.1) ~ S.1, col='red', lwd = 2)
legend(-5, 105, c('YR = 1', 'YR = 0.9', 'YR = 0.8', 'YR = 0.5', 'YR = 1.1', 'YR = 1.2', 'YR = 1.5'), pch=19, col = c('black', 'gray20', 'gray40',  'gray60', 'gray20', 'gray40',  'gray60'), bty = 'n')
dev.off()

### do this for Elahi et al data
data1 = read.csv("./elahidata.csv")

library(plyr)
tab1 <- ddply(data1, .(subSiteID), summarize, min(date.no))
tab2 <- ddply(data1, .(subSiteID), summarize, min(initialDate.no))
tab3 <- ddply(data1, .(subSiteID), summarize, max(date.no))
names(tab3) <- c('subSiteID', 'endDate.no')
data2 <- merge(data1, tab3, by = 'subSiteID')
data2$endRich <- ifelse(data2$date.no == data2$endDate.no, 'E','C')

data3 <- data2[data2$endRich == 'E',]
head(data3)
tab4 <- ddply(data3, .(subSiteID), summarize, max(rich))
names(tab4) <- c('subSiteID', 'endRich')
data4 <- merge(data2, tab4, by = 'subSiteID')

data4$dSR <- data4$endRich.y - data4$initialRich
hist(data4[data4$Trophic == '1',]$dSR, breaks = 40)
hist(data4[data4$Trophic == '2',]$dSR, breaks = 40)
hist(data4[data4$Trophic == '1',]$initialRich, breaks = 40)
hist(data4[data4$Trophic == '2',]$initialRich, breaks = 40)
plot(data4$dSR ~ data4$initialRich)

## for primary producers
b <- 0.26 # scaling exponent
S.1 <- seq(1, 150, 1) #species richness before loss
#S.2 number of species lost
#Y.R proportion biomass remaining after species loss
S.2 = function(Y, S.1) exp((1/b) * ( log(Y) + b*log(S.1) ))

pdf(file = "figure RBE alg.pdf", width = 6, height = 6)
plot(data4[data4$Trophic == '1',]$endRich.y ~ data4[data4$Trophic == '1',]$initialRich, pch = 19, ylim = c(0,40), xlim = c(0,40), xlab = 'Species Richness before change', ylab = 'Species Richness after change', main = ('Elahi data algae'))
abline(a = 1, b = 1, lwd = 2)
lines(S.2.func(0.9, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(0.8, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(0.5, S.1) ~ S.1, col='red', lwd = 2)
lines(S.2.func(1.1, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(1.2, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(1.5, S.1) ~ S.1, col='red', lwd = 2)
legend(-5, 105, c('YR = 1', 'YR = 0.9', 'YR = 0.8', 'YR = 0.5', 'YR = 1.1', 'YR = 1.2', 'YR = 1.5'), pch=19, col = c('black', 'gray20', 'gray40',  'gray60', 'gray20', 'gray40',  'gray60'), bty = 'n')
dev.off()

pdf(file = "figure deltaY Elahi alg.pdf", width = 4, height = 4)
hist(Y.func(data4[data4$Trophic == '1',]$endRich.y, data4[data4$Trophic == '1',]$initialRich), col = 'gray40',  breaks = 40, main = 'Change in function predicted for Elahi alg', xlab = 'Expected Biomass Change (Y.2/Y.1)', las = 1, cex = 1.4, xlim = c(0.5, 1.5))
abline(v = 0.9, lwd = 2)
abline(v = 1.1, lwd = 2)
dev.off()

## for grazers
b <- 0.47 # scaling exponent
S.1 <- seq(1, 100, 1) #species richness before loss
#S.2 number of species lost
#Y.R proportion biomass remaining after species loss
S.2 = function(Y, S.1) exp((1/b) * ( log(Y) + b*log(S.1) ))


pdf(file = "figure RBE herbs.pdf", width = 6, height = 6)
plot(data4[data4$Trophic == '2',]$endRich.y ~ data4[data4$Trophic == '2',]$initialRich, pch = 19, ylim = c(0,100), xlim = c(0,100), xlab = 'Species Richness before change', ylab = 'Species Richness after change', main = ('Elahi data herbivores'))
abline(a = 1, b = 1, lwd = 2)
lines(S.2.func(0.9, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(0.8, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(0.5, S.1) ~ S.1, col='red', lwd = 2)
lines(S.2.func(1.1, S.1) ~ S.1, col='gray20', lwd = 2)
lines(S.2.func(1.2, S.1) ~ S.1, col='dark blue', lwd = 2)
lines(S.2.func(1.5, S.1) ~ S.1, col='red', lwd = 2)
legend(-5, 105, c('YR = 1', 'YR = 0.9', 'YR = 0.8', 'YR = 0.5', 'YR = 1.1', 'YR = 1.2', 'YR = 1.5'), pch=19, col = c('black', 'gray20', 'gray40',  'gray60', 'gray20', 'gray40',  'gray60'), bty = 'n')
dev.off()

pdf(file = "figure deltaY Elahi herbs.pdf", width = 4, height = 4)
hist(Y.func(data4[data4$Trophic == '2',]$endRich.y, data4[data4$Trophic == '2',]$initialRich), col = 'gray40', breaks = 40, main = 'Change in function predicted for Elahi Herbivores', xlab = 'predicted proportional change in function', las = 1, cex = 1.4)
abline(v = 0.9, lwd = 2)
abline(v = 1.1, lwd = 2)
dev.off()



## plot of all scaling functions

pdf(file = "figure 1C.pdf", width = 6, height = 4)
plot(5*S.1^b ~ S.1, pch = '', xlim = c(1, 150), ylim = c(0, 80), xlab = 'Species richness (S)', ylab = 'estimated biomass (Y)')
b <- 0.26
lines(5*S.1^b ~ S.1, lwd = 2, col = 1)

b <- 0.47
lines(5*S.1^b ~ S.1, lwd = 2, col = 'gray40')

b <- 0.52
lines(5*S.1^b ~ S.1, lwd = 2, col = 'gray60')

legend(0, 80, c('Plants and Algae', 'Aq. Herbivores', 'Aq. Detritivores'), pch=19, col = c('black', 'gray40', 'gray60'), bty = 'n')
dev.off()
