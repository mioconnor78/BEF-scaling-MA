#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure 1
### Jan 20 2015: revised to explore relationship between Smax and random intercept coefficient in FM 
### Author: Mary O'Connor
#######################################################################


##### Producing experiment-level slope estimates #########
##########################################################
data <- SST4

rand.cat <- data.frame(cbind(as.numeric(as.character(data$Entry)), as.numeric(as.character(data$Study))))
names(rand.cat) <- c('Entry', 'Study')
rand.cat1 <- ddply(rand.cat, .(Entry), summarize, mean(Study)) 

rand.cat <- ddply(data, .(Entry, Study, Sys1, TG1, HigherT, restrt, unit.types2, logMaxTime, logSmaxc, LnTscale), summarize, mean(logY.rs))
#rand.cat$Syst <- rand.cat$Sys1
names(rand.cat) <- c('Mno', 'Study', 'Syst','TG1', 'HT', 'restrt','unit.types2', 'MaxTime', 'Smax', 'Tscale','meanlogY')
Entry.coefs <- data.frame(coef(modFM)$Entry)
Entry.coefs$Entry <- rownames(Entry.coefs)
#names(Mno.coefs)<-c('Int', 'LogS', 'System', 'HT', 'LogS.Sys1', 'LogS.HT', 'Mno1')
S <- cbind(rand.cat, Entry.coefs)

S$Sys.term <- ifelse(S$Syst == 'T', S$logSc.Sys1T, 0)
S$TG.term<-ifelse(S$TG1 == '2', S$logSc.TG12, 0)
S$TG.term<-ifelse(S$TG1 == '3', S$logSc.TG13, S$TG.term)
S$TG.term<-ifelse(S$TG1 == '4', S$logSc.TG14, S$TG.term)
S$units.term <- ifelse(S$unit.types2 == 'density', S$logSc.unit.types2density, 0)
S$units.term <- ifelse(S$unit.types2 == 'perc.cover', S$logSc.unit.types2perc.cover, S$units)
S$HT.term <- ifelse(S$HT == 'Y', S$logSc.HigherTY, 0)
S$Res.term <- ifelse(S$restrt == 'incr', S$logSc.restrtincr, 0)
S$Res.term <- ifelse(S$restrt == 'red', S$logSc.restrtred, S$Res.term)

#now add Study ranefs
St.ranefs <- data.frame(coef(modFM)$Study)
St.ranefs$Study <- rownames(St.ranefs)
St.ranefs1 <- data.frame(St.ranefs$Study, St.ranefs$logSc) # St.ranefs$X.Intercept.
#names(St.ranefs1) <- c('Study', 'logSc', 'Intercept')
S2 <- merge(S, St.ranefs1, by.x = 'Study', by.y = 'St.ranefs.Study', all= FALSE)
S <- S2
b <- as.numeric(fixef(modFM)[2])
S$slope <- S$logSc + S$Sys.term + S$TG.term + S$HT.term + S$units.term + S$Res.term + S$logSc.log.Tscale. + S$logSc.logSmaxc + S$logSc.log.MaxTscale...1. + (S$St.ranefs.logSc - b)
par(mar = c(4.5,4.5,3,2))
hist(S$slope, breaks = 40, col = 'gray', freq = TRUE, main = '', xlab = 'Estimated scaling coefficients (b)', xlim = c(-0.4, 1.4), ylim = c(0, 100), cex.lab = 1.2, axes = FALSE, ylab = 'Number of experiments') #
axis(1, at =c(-0.4, -0.2, 0, 0.2, 0.4, 0.6,0.8,1.0,1.2,1.4), lwd = 2, pos = 0)
axis(2, lwd = 2, pos = -0.4, las = 2)

summary(S$slope)  # ok, this is genreally working, though something is not quite right. the mean from this is 0.26, but the mean slope from the base model is 0.23. OK, the difference between these two is that the fixef for logS was counted twice - once in Mno and once in Study. that is now corrected by subtracting b in the final term above.
abline(v = mean(S$slope), lwd = 2)

### calculating CIs for the estimates:
n<-length(S$slope)
est<- mean(S$slope)
se <- sd((S$slope/sqrt(n)))
in.95 <- est + qt(c(0.025, 0.975), n-1)*se
abline(v = in.95[1], lwd = 2, lty = 2)
abline(v = in.95[2], lwd = 2, lty = 2)



# exploring correlation beteen Smax and ranefs
St.ranefs1

plot(S$X.Intercept. ~ S$Smax)
summary(lm(S$X.Intercept. ~ log(S$Smax)))

plot(S$X.Intercept. ~ as.numeric(S$Tscale))
summary(lm(S$X.Intercept. ~ as.numeric(S$Tscale)))

plot(S$X.Intercept. ~ as.numeric(S$Tscale))
summary(lm(S$X.Intercept. ~ as.numeric(S$MaxTime)))
