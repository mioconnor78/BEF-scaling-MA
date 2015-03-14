#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure 1
### Jan 06 2015; Author: Mary O'Connor
#######################################################################


##### Producing experiment-level slope estimates #########
##########################################################
data <- SST5
mod <- modBasic

rand.cat <- data.frame(cbind(as.numeric(as.character(data$Entry)), as.numeric(as.character(data$Study))))
names(rand.cat) <- c('Entry', 'Study')
rand.cat1 <- ddply(rand.cat, .(Entry), summarize, mean(Study)) 

### for modBasic
rand.cat <- ddply(data, .(Entry, Study), summarize, mean(logY.rs))
names(rand.cat) <- c('Mno', 'Study', 'meanlogY')
Entry.coefs <- data.frame(coef(mod)$Entry)
Entry.coefs$Entry <- rownames(Entry.coefs)
S <- cbind(rand.cat, Entry.coefs)

### for modBtrophic; skip to line 42 if using modBasic
mod <- modBtrophic
rand.cat <- ddply(data, .(Entry, Study, Sys1, TG1, HigherT, restrt), summarize, mean(logY.rs))
names(rand.cat) <- c('Mno', 'Study', 'Syst','TG1', 'HT', 'restrt', 'meanlogY')
Entry.coefs <- data.frame(coef(modBtrophic)$Entry)
Entry.coefs$Entry <- rownames(Entry.coefs)
S <- cbind(rand.cat, Entry.coefs)

S$Sys.term <- ifelse(S$Syst == '1', S$logSc.Sys1, 0)
S$TG.term<-ifelse(S$TG1 == '2', S$logSc.TG12, 0)
S$TG.term<-ifelse(S$TG1 == '3', S$logSc.TG13, S$TG.term)
S$TG.term<-ifelse(S$TG1 == '4', S$logSc.TG14, S$TG.term)
#S$units <- ifelse(S$Units == 'density', S$logS.unit.types2density, 0)
#S$units <- ifelse(S$Units == 'perc.cover', S$logS.unit.types2perc.cover, S$units)
S$HT.term <- ifelse(S$HigherT == 'Y', S$logSc.HT, 0)
#S$Res.term <- ifelse(S$restrt == 'incr', S$logS.restrtincr, 0)
#S$Res.term <- ifelse(S$restrt == 'red', S$logS.restrtred, S$Res.term)

#now add Study ranefs
St.ranefs <- data.frame(coef(mod)$Study)
St.ranefs$Study <- rownames(St.ranefs)
St.ranefs1 <- data.frame(St.ranefs$Study, St.ranefs$logSc)
S2 <- merge(S, St.ranefs1, by.x = 'Study', by.y = 'St.ranefs.Study', all= FALSE)
S <- S2
b <- as.numeric(fixef(mod)[2])

### for modBasic only:
S$slope <- S$logSc + S$logSc.log.Tscale. + (S$St.ranefs.logSc - b)

### for modBtrophic only: 
S$slope <- S$logSc + S$Sys.term + S$TG.term + S$HT.term + S$logSc.log.Tscale. + (S$St.ranefs.logSc - b)


## plot histogram of estimated slopes
par(mar = c(4.5,4.5,3,2))
hist(S$slope, breaks = 40, col = 'gray', freq = TRUE, main = '', xlab = 'Estimated scaling coefficients (b)', xlim = c(-0.4, 1.4), ylim = c(0, 100), cex.lab = 1.2, axes = FALSE, ylab = 'Number of experiments') #
axis(1, at =c(-0.4, -0.2, 0, 0.2, 0.4, 0.6,0.8,1.0,1.2,1.4), lwd = 2, pos = 0)
axis(2, lwd = 2, pos = -0.4, las = 2)

## identify extreme random effects (for comparison with caterpillar plots)
S$Study <- as.numeric(as.character(S$Study))
S$ranefs <- (S$St.ranefs.logSc - b) + (S$logSc - b)
hist(S$ranefs)
hist(S$St.ranefs.logSc - b)
hist(S$logSc - b)
S[(which(S$ranefs< -0.18)),1] # gives studies with very low ranefs
S[(which(S$ranefs> 0.22)),1] # gives studies with very low ranefs

names(S)
hist(S$X.Intercept.)
S$Int.ranef<-(S$St.ranefs)

mod <- modBtrophic
class(ranef(mod))
Entry.ranefs <- as.data.frame(ranef(mod)[1]) ## lmer::ranef does this too
Study.ranefs <- as.data.frame(ranef(mod)[2])
Entry.ranefs$Entry <- rownames(Entry.ranefs)
names(Entry.ranefs) <- c('Ent.rint', 'En.rslope', 'Entry')
head(Entry.ranefs)
Study.ranefs$Study <- rownames(Study.ranefs)
names(Study.ranefs) <- c('St.rint', 'St.slope', 'Study')
head(Study.ranefs)
S3 <- merge(S, Study.ranefs, by.x = 'Study', by.y = 'Study', all= FALSE)
S4 <- merge(S3, Entry.ranefs, by.x = 'Mno', by.y = 'Entry', all= FALSE)
S4$r.int <- S4$Ent.rint + S4$St.rint
S4$r.sl <- S4$St.slope + S4$En.rslope
hist(S4$r.int)
hist(S4$r.sl)
S4[(which(S4$r.int< -3.5)),1:2]


summary(S$slope)  # ok, this is genreally working, though something is not quite right. the mean from this is 0.26, but the mean slope from the base model is 0.23. OK, the difference between these two is that the fixef for logS was counted twice - once in Mno and once in Study. that is now corrected by subtracting b in the final term above.
abline(v = mean(S$slope), lwd = 2)

### calculating CIs for the estimates:
n<-length(S$slope)
est<- mean(S$slope)
se <- sd((S$slope/sqrt(n)))
in.95 <- est + qt(c(0.025, 0.975), n-1)*se
abline(v = in.95[1], lwd = 2, lty = 2)
abline(v = in.95[2], lwd = 2, lty = 2)
abline(v = est, lwd = 2, lty = 1)

x <- seq(0, 10, 1)
modx <- function(x)(yield~0.28*x)
plot9modx