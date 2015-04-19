#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### data formatting in preparation for mixed effects analysis
### we're going with this version; it processes the file with an ExptA column
### Mar 2015; Author: Mary O'Connor
#######################################################################

## Analysis of the relationship between SST and S for Cardinale et al BEF database
## Analysis of data without considering Taxon as an intercept predictor

## load data and libraries
library(plyr)
library(ggplot2)
library(reshape2)

### adapting lefcheck's code for data processing:
## I had to clean the metamaster file:
# removing extra levels in FinalT column and HigherT and consolidating to Y or N
# added HigherT values for Bret-Harte and Erskine, Gtim of 10 days to fungi studies
# converted FTG to numerical such that: (C = '3', P = '1', H = '2', D = '4', M = '5', O = '6'))

#Read in meta master
metamaster=read.csv("./BEF_MetaMaster_2011_08_29_exptA.csv")

#Subset data
metamaster2=ddply(metamaster,1,.progress="text",function(x) { 
  y=melt(x,id.vars=c(2:3,5,7:9,11,15,26,28,34:39,41),measure.vars=c(100:127)) 
  y$value=as.numeric(as.character(y$value))
  z=cbind(y[,1:17],richness=as.numeric(gsub("\\D","",y$variable)),value=y$value) 
  z=z[!is.na(z[,7]),] } ) 

# add mean vals for standardizing; Entry is the right group for this
mean.vals <- ddply(metamaster2, .(Entry), summarise, mean(value, na.rm = TRUE))
names(mean.vals) <- c('Entry', 'Mean.value')
metamaster.means <- merge(metamaster2,mean.vals, by.x = 'Entry', by.y = 'Entry')
metamaster.means$value.st <- metamaster.means$value/metamaster.means$Mean.value
head(metamaster.means)

#bring in restrt col from mary's file
mo<-read.csv("./input.HMM.stackn0unit23.csv", sep=",",header=T, na.strings="NA", fill=TRUE);
restrt <- ddply(mo, .(Entry, Mno, restrt, Study), summarize, mean(YEmono))
metamaster3 <- merge(restrt, metamaster.means, by.x = "Entry", by.y = "Entry", all = TRUE)

SST<-subset(metamaster3, metamaster.means$Ygen=='SST', select=1:25, drop=TRUE)
SST<-subset(SST, SST$value.st!='NaN', select=1:25, drop=TRUE) 
#SSTa<-subset(SST, SST$value!='NA', select=1:25, drop=TRUE) 

## simplifying Yunits: try to categorize units (e.g., n.density, m.density )
unit.type<-c('normalized mass', 'mass.normalized flux', 'mass.normalized flux', 'mass.normalized flux','mass.normalized flux', 'vol flux','mass flux', 'cover','rate','mass rate', 'normalized mass', 'mass.normalized flux','mass','mass', 'mass.normalized flux', 'vol flux', 'vol flux', 'cover', 'mass.normalized flux', 'mass.normalized flux', 'mass.normalized flux', 'normalized mass', 'normalized mass', 'mass.normalized flux','mass', 'mass.normalized flux','mass.normalized flux', 'normalized mass', 'mass flux', 'rate','rate','rate','mass.vol', 'mass', 'density', 'density', 'density', 'density', 'density','proportional change', 'mass.normalized flux','rate')
Y.units<-c(levels(SST$Yunits))
unit.types<-as.data.frame(cbind(Y.units, unit.type))
step1<-unit.types$unit.type[match(SST$Yunits, unit.types$Y.units)]
SST$unit.types<-as.factor(step1)

unit.type2<-c('biomass', 'flux', 'flux', 'flux','flux', 'flux','flux', 'perc.cover','flux','flux', 'biomass', 'flux','biomass','biomass', 'flux','flux', 'flux', 'biomass', 'flux', 'flux', 'flux', 'biomass', 'biomass', 'flux', 'biomass', 'flux','flux', 'biomass','flux','flux','flux', 'biomass', 'biomass','biomass','density', 'density', 'density', 'density', 'density', 'proportional change', 'flux', 'flux')
unit.types2<-as.data.frame(cbind(Y.units, unit.type2))
step2<-unit.types2$unit.type2[match(SST$Yunits, unit.types2$Y.units)]
SST$unit.types2<-as.factor(step2)

## removing levels or values not relevant for this analysis
SST<-subset(SST, SST$TDBU=='TD', select=1:27, drop=TRUE) # removing 'bottom up' studies
SST1<-SST[-which(SST$Entry=='616'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations
SST1<-SST1[-which(SST1$Entry=='617'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations
SST1<-SST1[-which(SST1$Entry=='618'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations
SST1<-SST1[-which(SST1$Entry=='619'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations 
SST1<-SST1[-which(SST1$Entry=='250'),] # removing entry for Mikola 1998 because a) I can't understand where it came from when i read the paper, and b) the data from the figures in the paper is present in other entries
#SST1<-SST1[-which(SST1$Entry=='246'),] # Mikola 1998, ditto entry 250; this is gone now from removing NaNs above.
#SST2<-subset(SST2, SST2$value!= '0', select=1:27, drop=TRUE) 
SST2<-subset(SST1, SST1$Slevels>1, select=1:27, drop=TRUE) 
SST2<-subset(SST2, SST2$HigherT!="", select=1:27, drop=TRUE)
SST2<-subset(SST2, SST2$HigherT!=".", select=1:27, drop=TRUE)
# get rid of Tscale vals = 0
SST2<-subset(SST2, SST2$Tscale!="", select=1:27, drop=TRUE) 
SST2<-subset(SST2, SST2$Yunits!='proportional change', select=1:27, drop=TRUE) 
SST2 <- SST2[which(SST2$Yunits!='rate'),]


### need to estimate the duration of each experiment:
maxTime <- ddply(SST2, .(Entry, Tscale), summarise, TFinal = 'Y')
names(maxTime) <- c('Entry', 'MaxTscale', 'Finaltime')
merge(SST2, maxTime, by.x = "Entry", by.y = "Entry", all.x = TRUE, all.y = FALSE) -> merged
dim(merged)
names(merged)
head(merged)
SST2<-merged
SST2$logMaxTime <- log(as.numeric(as.character(SST2$MaxTscale)))

## transform columns
SST2$logY <- log(SST2$value)
SST2$logS <- log(SST2$richness)
SST2$TG1 <- as.factor(SST2$FTG)
SST2$Tscale <- as.numeric(as.character(SST2$Tscale))
SST2$Smax <- as.numeric(as.character(SST2$Smax))
SST2$MaxTscale <- as.numeric(as.character(SST2$MaxTscale))
SST2$Entry <- as.factor((SST2$Entry))
SST2$Study <- as.factor((SST2$Study))
SST2 <- SST2[,-(4)]

## rescaling entries if their values are extremely low or high
try1 <- ddply(SST2, .(Entry, Study, Yunits), summarise, min(value))
try2 <- ddply(SST2, .(Entry), summarise, max(value))
names(try1) <- c('Entry', 'Study','units','minval')
names(try2) <- c('Entry', 'maxval')
merge(try1, try2, by.x = "Entry", by.y = "Entry") -> try3
try3$convert.min <- ifelse(try3$maxval < 10, '1', '0')
try3$convert.max <- ifelse(try3$minval > 10000, '1', '0')
sorted <- try3[(order(try3$minval)),]
try3<-try3[,-(2:3)]
merge(SST2, try3, by.x = "Entry", by.y = "Entry") -> SST4

SST4$values.rs <- as.numeric(as.character(ifelse(SST4$convert.min == 1, SST4$value*1000, SST4$value)))
SST4$values.rs <- as.numeric(as.character(ifelse(SST4$convert.max == 1, SST4$value/1000, SST4$values.rs)))
SST4$logY.rs <- log(SST4$values.rs)
head(SST4)
plot(SST4$logY.rs ~ SST4$logS, main = 'SST4.rs2')

## remove outliers based on previous analysis using visual inspection of plot(modBasic)
dim(SST4)
SST4 <- subset(SST4, SST4$Study!=83, select=1:38, drop=TRUE)
SST4 <- subset(SST4, SST4$Study!=177, select=1:38, drop=TRUE) 
#SST4 <- subset(SST4, SST4$Mno!=796, select=1:38, drop=TRUE) # based on looking at residuals of individual regressions, this one is an extreme outlier (below)
#SST4 <- subset(SST4, SST4$Mno!=826, select=1:38, drop=TRUE) # searching for the outlier in plot(modF1)


SST4$logYst <- log(SST4$value.st)
## center the regressor
SST4$logSc <- SST4$logS - log(8)
SST4$lnTc <- log(SST4$Tscale) - log(mean(SST4$Tscale))
  #log(mean(SST4$richness))
  #mean(SST4$logS)
SST4$logSmax <- log(SST4$Smax)
SST4$logSmaxc <- SST4$logSmax - median(SST4$logSmax)

plot(SST4$logY.rs ~ SST4$logSc, main = 'SST4.rs2')
## the units column will be wrong for rescaled values, but in the model we use 'unit.types', and that class should still be fine.
## upon inspection, I can see that some studies (e.g., 8) will have some rescaled values and some not rescaled, which would bring the intercepts together. shouldn't be a problem.

#remove carnivores
SST5 <- subset(SST4, SST4$TG1!="3", select=1:40, drop=TRUE) 

########## INITIAL DATA PREP COMPLETE #########

## removing studies with extreme random effects to see if we can get rid of that coefficient:
SST5 <- subset(SST5, SST5$Study!=147, select=1:40, drop=TRUE) 
SST5 <- subset(SST5, SST5$Study!=87, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=174, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=66, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=103, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=121, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=127, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=137, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=168, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=183, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=30, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=155, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=80, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=157, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=180, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=192, select=1:40, drop=TRUE)
SST5 <- subset(SST5, SST5$Study!=84, select=1:40, drop=TRUE)
## by this time, the correlation between study and entry is going up, but the range of entry random effects is declined (either study 155 or 80 had really high intercept ranefs), and the model coefs are shifting for modBtrophic. 


### DATA PROCESSING IS COMPLETE ###

## testing for correlations among model parameters:
logSc*log(Tscale) + logSc*unit.types2 + logSc*logSmaxc + logSc*log(MaxTscale+1)

#plot(log(as.numeric(SST4$Tscale)) ~ SST4$logSmaxc, main = 'logTscale ~ Smaxc')
#summary(lm(log(as.numeric(SST4$Tscale)) ~ SST4$logSmaxc))

plot(log(as.numeric(log(SST4$MaxTscale+1))) ~ SST4$logSmaxc, main = 'logTscale ~ Smaxc')
summary(lm(log(as.numeric(SST4$MaxTscale+1)) ~ SST4$logSmaxc))

plot(SST4$logYst~ log(as.numeric(SST4$Tscale)), main = 'logYst ~ Tscale')
summary(lm(SST4$logYst~ log(as.numeric(SST4$Tscale))))

plot(SST4$logYst~ log(as.numeric(log(SST4$MaxTscale+1))), main = 'logYst ~ MaxTscale')
summary(lm(SST4$logYst~ log(as.numeric(log(SST4$MaxTscale+1)))))

plot(SST4$logYst~ SST4$logSmaxc, main = 'logYst ~ Smaxc')
summary(lm(SST4$logYst~ SST4$logSmaxc))


## some exploration
plot(SST4[(SST4$TG1=='3'),]$logY.rs ~ SST4[(SST4$TG1=='3'),]$logS)
plot(SST4$value.st ~ SST4$logS, main = 'SST4.rs2')
plot(SST4$logYst ~ SST4$logS, main = 'SST4.rs2')
plot(SST4$logYst ~ SST4$logSc, main = 'logY.st ~ logSc')
plot(SST4$logYst ~ SST4$logSmax, main = 'logY.st ~ lnSmax')

hist(SST4$Smax)
hist(SST4$logSmax)
hist(SST4$value.st)
hist(SST4$logYst)
hist(SST4$lnTscale)
#hist(SST4[(SST4$logYst>-2),]$logYst)

plot(SST4$logYst ~ SST4$logSmax, main = 'logY.st ~ Smax')
summary(lm(SST4$logYst ~ SST4$logSmax))

SST4[(SST4$value.st <= .3),]




## data summary for Table S1
length(unique(SST5$Entry))

ddply(SST5, .(Sys1), summarize, length(unique(Entry)))
ddply(SST5, .(Sys1), summarize, length(unique(ExptA)))
ddply(SST5, .(Sys1), summarize, length(unique(Study)))

ddply(SST5, .(TG1), summarize, length(unique(Entry)))
ddply(SST5, .(TG1), summarize, length(unique(ExptA)))
ddply(SST5, .(TG1), summarize, length(unique(Study)))

ddply(SST5, .(HigherT), summarize, length(unique(Entry)))
ddply(SST5, .(HigherT), summarize, length(unique(ExptA)))
ddply(SST5, .(HigherT), summarize, length(unique(Study)))

ddply(SST5, .(unit.types2), summarize, length(unique(Entry)))
ddply(SST5, .(unit.types2), summarize, length(unique(ExptA)))
ddply(SST5, .(unit.types2), summarize, length(unique(Study)))

ddply(SST5, .(restrt), summarize, length(unique(Entry)))
ddply(SST5, .(restrt), summarize, length(unique(ExptA)))
ddply(SST5, .(restrt), summarize, length(unique(Study)))

res.test <- ddply(restrt, .(Study), summarize, length(unique(restrt)))
res.test[(res.test$..1 > 1),]

res.test <- ddply(SST5, .(ExptA, Ref), summarize, (unique(restrt)))


refs <- ddply(SST, .(Ref), summarize, length(Entry))
write.csv(refs, 'refs.csv')


## exploring error sturcture

plot(SST4$values.rs ~ SST4$richness)
means <- ddply(SST4, .(Entry), summarise, mean(values.rs))
names(means) <- c('Entry', 'means')
vars <- ddply(SST4, .(Entry), summarise, var(values.rs))
names(vars) <- c('Entry', 'vars')
means.vars <- merge(means, vars, by.x = 'Entry', by.y = 'Entry')
par(mfrow = c(2,2))
plot(means.vars$means, means.vars$vars)
plot(means.vars$means, means.vars$vars, xlim = c(0, 2000), ylim = c(0, 4000))
plot(means.vars$means, means.vars$vars, xlim = c(0, 1000), ylim = c(0, 4000))
lines(lowess(means.vars$means), col = 2)
plot(mean(SST4$values.rs ~ SST4$richness)
     
rich.range <- ddply(SST4, .(Entry), summarise, max(richness))


x <- 0:100

y <- rnorm(101, x^0.2, .1)

plot(x,y)

summary(lm(log(y) ~ x))
summary(lm(log(y/mean(y)) ~ x))
