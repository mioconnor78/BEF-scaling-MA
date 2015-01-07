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

write.csv(SST2, 'SST2new.csv')

## troubleshooting linear models:
study.test <- ddply(SST2, .(Study, Ref), summarise, (count(Entry)))
study.test2 <- ddply(SST2, .(Study, Entry), summarise, (mean(Entry)))

modBasic <- lm(logY ~ logS, data=data, na.action=na.omit)

SST2[SST2$logY=='Inf',]
SST2[SST2$logY=='-Inf',]


## 1.1 load and format the datafile
mo<-read.csv("/Users/maryo/Dropbox/nceas_bdef/Scaling Relationship/metaanalysis of B/analysis HMM 2012/input.HMM.stackn0unit23.csv", sep=",",header=T, na.strings="NA", fill=TRUE);
mo1<-mo[,-(2:49)]
head(mo1)
names(mo1)<-c('Mno',seq(1,72,1))
stacked<-stack(mo1)
length(mo1[,2])
stacked<-stacked[-(1:1415),]
### Mno is a unique identifer for each row in the meta-master
Mno<-mo1$Mno
stacked$Mno<-Mno
names(stacked)<-c('Y','richness', 'Mno')
stacked$richness.n<-rep(1:72, each=length(stacked[,1])/72)

mo2<-mo[,-(49:121)]
both<-merge(stacked, mo2, by.x='Mno', by.y = 'Mno', all=T) 

## 1.2 cleaning up
both$Y<-as.numeric(as.character(both$Y))
both$richness.n<-as.numeric(as.character(both$richness.n)) 
both$logS<-log(both$richness.n)
both$Smax<-as.factor(both$Smax)
both$Mno<-as.numeric(as.character(both$Mno))
both$Smax1<-as.numeric(as.character(both$Smax))
both$Sys1<-as.numeric(as.character(both$Sys1))
both$TG1<-as.factor((both$TG))

## 1.3 Removing row types that are problematic for analysis
## i'm going to remove proportional change responses; they are often negative, which presents a problem for log transformation.
both<-subset(both, both$HigherT!="", select=1:54, drop=TRUE)
## We can't really inlcude spatial scale or body size; we loose all the non-plant trophic levels. so options here would be to go get the body size data, or just do that analysis for plants (though i suspect it's been done)
# get rid of Tscale vals = 0
both<-subset(both, both$Tscale!="", select=1:54, drop=TRUE) 
both<-subset(both, both$Yunits!='proportional change', select=1:54, drop=TRUE) 
both<-subset(both, both$Y!='NA', select=1:54, drop=TRUE) 
negs<-both[which(both$Y<=0),]
## only 3 negatives and 23 Y=0, and they are not SST vals, so i'm going to remove them for now
both<-subset(both, both$Y>0, select=1:54, drop=TRUE) 
both$logY<-log(both$Y) 

## 1.5 make datafile of just SST and top down only
SSTa<-subset(both, both$Ygen=='SST', select=1:57, drop=TRUE)
#SST<-subset(SST, SST$FinalT=='Y', select=1:57, drop=TRUE) # this gives us a dataset using only final years.
SSTa<-subset(SSTa, SSTa$TDBU=='TD', select=1:57, drop=TRUE) # removing 'bottom up' studies
SSTa<-subset(SSTa, SSTa$Slevels>1, select=1:57, drop=TRUE) 

#restrt <- ddply(SSTa, .(Entry, Mno, restrt), summarize, mean(logY))
#SSTb <- merge(restrt, SST2, by.x = "Entry", by.y = "Entry", all = TRUE)

# explored these as outliers, but they are not, or I ended up removing the whole study:
#SST2<-subset(SST, SST$Study!=27, select=1:57, drop=TRUE) # the study that contains 796 and other extreme values
#SST2<-subset(SST2, SST2$Mno!=543, select=1:57, drop=TRUE) # removing the outlying Mno from study 177; nope, both are outliers.
#SST2 <- subset(SST2, SST2$Mno!=962, select=1:57, drop=TRUE) # this is the offender from study 83

#re: outliers: could remove both studies 177 and 83; OR, keep both and remove just the lowest value from study 83 (962 seems to be the offender)



SST2nc <- subset(SST2, SST2$FTG!="C", select=1:60, drop=TRUE) # the problematic detritivore study

