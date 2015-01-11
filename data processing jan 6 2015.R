#######################################################################
### Modeling uncertainty in the biodiversity ecosystem functioning relationship
### O'Connor, Gonzalez, et al
### data formatting in preparation for mixed effects analysis
### Jan 22 2014; Author: Mary O'Connor
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
metamaster=read.csv("./BEF_MetaMaster_2011_08_29.csv")

#Subset data
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
SST<-subset(SST, SST$TDBU=='TD', select=1:25, drop=TRUE) # removing 'bottom up' studies
SST1<-SST[-which(SST$Entry=='616'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations
SST1<-SST1[-which(SST1$Entry=='617'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations
SST1<-SST1[-which(SST1$Entry=='618'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations
SST1<-SST1[-which(SST1$Entry=='619'),] # removing douglass et al measurements of predator biomass for grazer diversity manipulations 
SST1<-SST1[-which(SST1$Entry=='250'),] # removing entry for Mikola 1998 because a) I can't understand where it came from when i read the paper, and b) the data from the figures in the paper is present in other entries
SST1<-SST1[-which(SST1$Entry=='246'),] # Mikola 1998, ditto entry 250
SST2<-subset(SST1, SST1$value!='NA', select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$value!= '0', select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$Slevels>1, select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$HigherT!="", select=1:25, drop=TRUE)
SST2<-subset(SST2, SST2$HigherT!=".", select=1:25, drop=TRUE)
# get rid of Tscale vals = 0
SST2<-subset(SST2, SST2$Tscale!="", select=1:25, drop=TRUE) 
SST2<-subset(SST2, SST2$Yunits!='proportional change', select=1:25, drop=TRUE) 
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
try3$convert.min <- ifelse(try3$maxval < 1, '1', '0')
try3$convert.max <- ifelse(try3$minval > 22000, '1', '0')
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
SST4 <- subset(SST4, SST4$Study!=177, select=1:36, drop=TRUE) 
SST4 <- subset(SST4, SST4$Mno!=796, select=1:36, drop=TRUE) # based on looking at residuals of individual regressions, this one is an extreme outlier (below)
SST4 <- subset(SST4, SST4$Mno!=826, select=1:36, drop=TRUE) # searching for the outlier in plot(modF1)
SST4 <- subset(SST4, SST4$Study!=83, select=1:36, drop=TRUE)

## the units column will be wrong for rescaled values, but in the model we use 'unit.types', and that class should still be fine.
## upon inspection, I can see that some studies (e.g., 8) will have some rescaled values and some not rescaled, which would bring the intercepts together. shouldn't be a problem.


#remove carnivores
SST5 <- subset(SST4, SST4$TG1!="3", select=1:36, drop=TRUE) 

## some exploration
plot(SST4[(SST4$TG1=='3'),]$logY.rs ~ SST4[(SST4$TG1=='3'),]$logS)


## data summary
length(unique(SST4$Entry))

ddply(SST4, .(Sys1), summarize, length(unique(Entry)))
ddply(SST4, .(Sys1), summarize, length(unique(Study)))
ddply(SST4, .(TG1), summarize, length(unique(Entry)))
ddply(SST4, .(TG1), summarize, length(unique(Study)))
ddply(SST4, .(HigherT), summarize, length(unique(Entry)))
ddply(SST4, .(HigherT), summarize, length(unique(Study)))
ddply(SST4, .(unit.types2), summarize, length(unique(Entry)))
ddply(SST4, .(unit.types2), summarize, length(unique(Study)))
ddply(SST4, .(restrt), summarize, length(unique(Entry)))
ddply(SST4, .(restrt), summarize, length(unique(Study)))


refs <- ddply(SST, .(Ref), summarize, length(Entry))
write.csv(refs, 'refs.csv')

