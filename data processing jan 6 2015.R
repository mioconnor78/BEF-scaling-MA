#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### data formatting in preparation for mixed effects analysis
### Jan 22 2014; Author: Mary O'Connor
#######################################################################

## Analysis of the relationship between SST and S for Cardinale et al BEF database
## Analysis of data without considering Taxon as an intercept predictor
## same as other Aug 24th, except in this version I will remove more outliers

## load data and libraries
library(plyr)
library(ggplot2)
library(reshape2)

### adapting lefcheck's code for data processing:
##  I had to clean the metamaster file:
# removing extra levels in FinalT column and HigherT
# added HigherT values for Bret-Harte and Erskine, Gtim of 10 days to fungi studies


#Read in meta master
metamaster=read.csv("./BEF_MetaMaster_2011_08_29.csv")

#Subset data and create response as proportional change in functioning with each increase in richness

metamaster2=ddply(metamaster,1,.progress="text",function(x) { 
  y=melt(x,id.vars=c(2:4,6:8,10,14,25,27,33:38,40),measure.vars=c(99:126)) 
  y$value=as.numeric(as.character(y$value))
  z=cbind(y[,1:17],richness=as.numeric(gsub("\\D","",y$variable)),value=y$value) 
  z=z[!is.na(z[,7]),] } ) 

#bring in restrt col from mary's file
mo<-read.csv("./input.HMM.stackn0unit23.csv", sep=",",header=T, na.strings="NA", fill=TRUE);
restrt <- ddply(mo, .(Entry, Mno, restrt), summarize, mean(YEmono))
metamaster3 <- merge(restrt, metamaster2, by.x = "Entry", by.y = "Entry", all = TRUE)

SST<-subset(metamaster3, metamaster2$Ygen=='SST', select=1:23, drop=TRUE)

## 1.4 simplifying Yunits: try to categorize units (e.g., n.density, m.density )
### this is amazing what i did here; i often want to do this but i forget that i already did it.
unit.type<-c('normalized mass', 'mass.normalized flux', 'mass.normalized flux', 'mass.normalized flux','mass.normalized flux', 'vol flux','mass flux', 'cover','rate','mass rate', 'normalized mass', 'mass.normalized flux','mass','mass', 'mass.normalized flux', 'vol flux', 'vol flux', 'cover', 'mass.normalized flux', 'mass.normalized flux', 'mass.normalized flux', 'normalized mass', 'normalized mass', 'mass.normalized flux','mass', 'mass.normalized flux','mass.normalized flux', 'normalized mass', 'mass flux', 'rate','rate','rate','mass.vol', 'mass', 'density', 'density', 'density', 'density', 'density','proportional change', 'mass.normalized flux','rate')
Y.units<-c(levels(SST$Yunits))
unit.types<-as.data.frame(cbind(Y.units, unit.type))
step1<-unit.types$unit.type[match(SST$Yunits, unit.types$Y.units)]
SST$unit.types<-as.factor(step1)

unit.type2<-c('biomass', 'flux', 'flux', 'flux','flux', 'flux','flux', 'perc.cover','flux','flux', 'biomass', 'flux','biomass','biomass', 'flux','biomass', 'biomass', 'flux', 'flux', 'perc.cover', 'flux', 'flux', 'flux', 'biomass', 'biomass', 'flux','biomass', 'flux','flux','flux','flux', 'biomass', 'biomass','biomass','density', 'density', 'density', 'density', 'density', 'proportional change', 'flux', 'flux')
unit.types2<-as.data.frame(cbind(Y.units, unit.type2))
step2<-unit.types2$unit.type2[match(SST$Yunits, unit.types2$Y.units)]
SST$unit.types2<-as.factor(step2)

SST<-subset(SST, SST$TDBU=='TD', select=1:25, drop=TRUE) # removing 'bottom up' studies
SST2<-subset(SST, SST$Slevels>1, select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$HigherT!="", select=1:25, drop=TRUE)
SST2<-subset(SST2, SST2$HigherT!=".", select=1:25, drop=TRUE)
# get rid of Tscale vals = 0
SST2<-subset(SST2, SST2$Tscale!="", select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$Yunits!='proportional change', select=1:25, drop=TRUE) 
SST2 <- SST[which(SST2$Yunits!='rate'),]
SST2<-subset(SST2, SST2$value!='NA', select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$value!= '0', select=1:25, drop=TRUE) 

### need to estimate the duration of each experiment:
maxTime <- ddply(SST2, .(Entry, Tscale), summarise, TFinal = 'Y')
names(maxTime) <- c('Entry', 'MaxTscale', 'Finaltime')
merge(SST2, maxTime, by.x = "Entry", by.y = "Entry", all.x = TRUE, all.y = FALSE) -> merged
dim(merged)
names(merged)
head(merged)
SST2<-merged
SST2$logMaxTime <- log(as.numeric(as.character(SST2$MaxTscale)))

SST2$logY <- log(SST2$value)
SST2$logS <- log(SST2$richness)
SST2$TG1 <- as.factor(SST2$FTG)
SST2$Tscale <- as.numeric(as.character(SST2$Tscale))
SST2$Smax <- as.numeric(as.character(SST2$Smax))
SST2$MaxTscale <- as.numeric(as.character(SST2$MaxTscale))
SST2$Entry <- as.integer((SST2$Entry))
SST2$Study <- as.factor((SST2$Study))

## outliers based on previous analysis
dim(SST2)
SST2 <- subset(SST2, SST2$Study!=177, select=1:31, drop=TRUE) # the outlier in modF1 results, determined by idenfitying the row of the residual value, in a dataframe
SST2 <- subset(SST2, SST2$Study!=83, select=1:31, drop=TRUE)
SST2 <- subset(SST2, SST2$Mno!=796, select=1:31, drop=TRUE) # based on looking at residuals of individual regressions, this one is an extreme outlier (below)
SST2 <- subset(SST2, SST2$Mno!=826, select=1:31, drop=TRUE) # searching for the outlier in plot(modF1)
SST2 <- SST2[,-(4)]

head(SST2)
write.csv(SST2, 'SST2new.csv')

## rescaling entries if their values are extremely low or high
try1 <- ddply(SST2, .(Entry), summarise, min(value))
try2 <- ddply(SST2, .(Entry), summarise, max(value))
names(try1) <- c('Entry', 'minval')
names(try2) <- c('Entry', 'maxval')
merge(try1, try2, by.x = "Entry", by.y = "Entry") -> try3
try3$convert.min <- ifelse(try1$minval < 1, '1', '0')
try3$convert.max <- ifelse(try2$maxval > 22000, '1', '0')
merge(SST3, try3, by.x = "Entry", by.y = "Entry") -> SST4
head(SST4)

# great. Now need to take all entries with a 1 and multiply their vals by 1000.
SST4$values.rs <- as.numeric(as.character(ifelse(SST4$convert == 1, SST4$value*1000, SST4$value)))
SST4$values.rs <- as.numeric(as.character(ifelse(SST4$convert.max == 1, SST4$value/1000, SST4$values.rs)))
SST4$logY.rs <- log(SST4$values.rs)
head(SST4)
plot(SST4$logY.rs ~ SST4$logS, main = 'SST4.rs2')

#extra stuff
SST2nc <- subset(SST2, SST2$FTG!="C", select=1:60, drop=TRUE) # the problematic detritivore study

