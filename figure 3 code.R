#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### figure 1
### Jan 06 2015; Author: Mary O'Connor
#######################################################################


##### Producing experiment-level slope estimates #########
##########################################################
data <- SST5

mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod <- mod4

rand.cat <- ddply(data, .(Entry, Study, ExptA, Sys1, TG1), summarize, mean(logY.rs))
names(rand.cat) <- c('Entry', 'Study', 'ExptA','Syst','TG1', 'meanlogY')
Entry.coefs <- data.frame(coef(mod)$Entry) #this gives us the logS for entry. so we'd need to add only the expt and study level later. let's try that.
Entry.coefs$Entry <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by = 'Entry', all = FALSE)

## constructing predicted slopes
S$Sys.term <- ifelse(S$Syst == 'T', S$logSc.Sys1T, 0)
S$TG.term <- ifelse(S$TG1 == '2', S$logSc.TG12, 0)
S$TG.term <- ifelse(S$TG1 == '4', S$logSc.TG14, S$TG.term)
S$TG.term <- ifelse((S$Syst=='T' & S$TG1=='4'), (S$logSc.TG14 + S$logSc.Sys1T.TG14), S$TG.term)

#now add Study and ExpA ranefs
St.ranefs       <- data.frame(coef(mod)$Study)
St.ranefs$Study <- rownames(St.ranefs)
Ex.ranefs       <- data.frame(coef(mod)$ExptA)
Ex.ranefs$ExptA <- rownames(Ex.ranefs)
#Ent.ranefs      <- data.frame(coef(mod)$Entry)
#Ent.ranefs$Entry <- rownames(Ent.ranefs)

St.ranefs1 <- data.frame(St.ranefs$Study, St.ranefs$logSc)
Ex.ranefs1 <- data.frame(Ex.ranefs$ExptA, Ex.ranefs$logSc)
#Ent.ranefs1 <- data.frame(Ent.ranefs$Entry, Ent.ranefs$logSc)

S2 <- merge(S, St.ranefs1, by.x = 'Study', by.y = 'St.ranefs.Study', all= FALSE)
S3 <- merge(S2, Ex.ranefs1, by.x = 'ExptA', by.y = 'Ex.ranefs.ExptA', all = FALSE)
#S3 <- merge(S3.1, Ent.ranefs1, by.x = 'Entry', by.y = 'Ent.ranefs.Entry', all = FALSE)
S <- S3
b <- as.numeric(fixef(mod)[2])



### for modBtrophic only: 
S$slope <- S$logSc + S$Sys.term + S$TG.term + (S$St.ranefs.logSc - b) + (S$Ex.ranefs.logSc - b) 
# + (S$Ent.ranefs.logSc - b)

summary(S$slope)

# checking some things:
length(S[S$Sys.term=='0',]$Sys.term)
length(S[S$TG.term!='0',]$TG.term)
plot(SST5[(SST5$Entry == '274'),]$logY.rs ~ SST5[SST5$Entry == '274',]$logSc, ylim = c(0, 10))
abline(4.2, S[(S$Entry == '274'),]$slope)

## plot histogram of estimated slopes
#pdf(file='Figure3.pdf', width = 5, height = 3)
par(mar = c(4.5,4.5,3,2))
hist(S$slope, breaks = 40, col = 'gray', freq = TRUE, main = '', xlab = 'Estimated scaling coefficients (b)', xlim = c(-0.4, 1.4), ylim = c(0, 100), cex.lab = 1.2, axes = FALSE, ylab = 'Number of entries') #
axis(1, at =c(-0.4, -0.2, 0, 0.2, 0.4, 0.6,0.8,1.0,1.2,1.4), lwd = 2, pos = 0)
axis(2, lwd = 2, pos = -0.4, las = 2)

### calculating CIs for the estimates:
n<-length(S$slope)
est<- mean(S$slope)
se <- sd((S$slope/sqrt(n)))
in.95 <- est + qt(c(0.025, 0.975), n-1)*se
abline(v = in.95[1], lwd = 2, lty = 2)
abline(v = in.95[2], lwd = 2, lty = 2)
abline(v = est, lwd = 2, lty = 1)

dev.off()

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


#this code used to be at the top to easily allow me to make this figure withi mod2F or mod3F

#rand.cat <- data.frame(cbind(as.numeric(as.character(data$Entry)), as.numeric(as.character(data$ExptA)), as.numeric(as.character(data$Study))))
#names(rand.cat) <- c('Entry', 'ExptA', 'Study')
#rand.cat1 <- ddply(rand.cat, .(Entry, ExptA), summarize, mean(Study)) 

### for modBasic
rand.cat <- ddply(data, .(Entry, ExptA, Study), summarize, mean(logY.rs))
names(rand.cat) <- c('Entry','ExptA', 'Study', 'meanlogY')
Entry.coefs <- data.frame(coef(mod)$Entry)
Entry.coefs$Entry <- rownames(Entry.coefs)

S <- cbind(rand.cat, Entry.coefs)

### for modBasic only:
S$slope <- S$logSc + S$logSc.log.Tscale. + (S$St.ranefs.logSc - b)


## this code could be added to the slope estimation section to include other model terms:

#S$units <- ifelse(S$Units == 'density', S$logS.unit.types2density, 0)
#S$units <- ifelse(S$Units == 'perc.cover', S$logS.unit.types2perc.cover, S$units)
#S$HT.term <- ifelse(S$HigherT == 'Y', S$logSc.HT, 0)
#S$Res.term <- ifelse(S$restrt == 'incr', S$logS.restrtincr, 0)
#S$Res.term <- ifelse(S$restrt == 'red', S$logS.restrtred, S$Res.term)