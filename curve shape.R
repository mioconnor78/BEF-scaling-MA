########################################################
##
## Code to analyze the Cardinale et al. BEF database
## and evalute curve shape for Y/Y(avg. monoculture)
## using a nonlinear mixed model apporach and AIC tables
##
## By Jon Lefcheck& Jarrett Byrnes
## modified Jan 06 2015 by Mary O'Connor to use data for scaling meta-analysis
##
########################################################

library(AICcmodavg)
library(nlme)
library(plyr)
library(reshape2)
library(ggplot2)

#Read in meta master
metamaster=read.csv("./BEF_MetaMaster_2011_08_29.csv")

#Subset data and create response as proportional change in functioning with each increase in richness

metamaster2=ddply(metamaster,1,.progress="text",function(x) { 
  y=melt(x,id.vars=c(2:4,7:8),measure.vars=c(99:126)) 
  y$value=as.numeric(as.character(y$value))
  z=cbind(y[,1:5],richness=as.numeric(gsub("\\D","",y$variable)),value=y$value/y[y$variable=="Y1","value"]) 
  z=z[!is.na(z[,7]),] } ) 



#Remove all studies with <2 levels of richness
metamaster.reduced=ddply(metamaster2,1:6,function(x) if(length(x$value)<3) NULL else x)



#Overall form

#Write function to fit models corresponding to different functional forms
modFit=function(df) {
  df=df[!is.na(df$value) & !is.infinite(df$value),]
  df=groupedData(value~richness|Ref,data=df)
  #Fit different models
  Null=nlme(value~a,fixed=a~1,random=~a~1,start=c(a=-1),control=nlmeControl(tolerance=1e-04),data=df)
  Linear=nlme(value~a+b*richness,fixed=a+b~1,random=~a+b~1,start=c(a=1.5,b=1),data=df)
  Logarithmic=nlme(value~a+b*log(richness),fixed=a+b~1,random=~a+b~1,start=c(a=1,b=1),data=df)
  Power=nlme(value~a*richness^b,fixed=a+b~1,random=~a+b~1,start=c(a=0.18,b=2),data=df)
  #Exponential=nlme(value~exp(a+b*richness),fixed=a+b~1,random=~a+b~1,start=c(a=1,b=1),data=df)
  Saturating=nlme(value~(max(value)*richness)/(k+richness),fixed=k~1,random=k~1,start=c(k=0.5),data=df)
  #Return models in list
  return(list(Null=Null,
              Linear=Linear, 
              Logarithmic=Logarithmic, 
              Power=Power, 
              #Exponential=Exponential, 
              Saturating=Saturating) ) }

#Write function to extract AIC values and weights from the model list obtained from function modFit
getAICtab=function(modList) {
  mixedAICs=data.frame(AIC=round(sapply(modList,AICc),2))
  mixedAICs=within(mixedAICs, {
    deltaAIC=round(AIC-min(AIC),1)
    modLik=exp(-0.5*deltaAIC)
    AICweight=round(modLik/sum(modLik),3) } ) } 

#Run for reduced dataset and subset by Ycat
mods=dlply(SST4,"Ygen",modFit)


#Show the AIC Table of Results
aicVals <- lapply(mods,getAICtab)
aicVals <- ldply(aicVals)
aicVals <- ldply(aicVals,function(i) { cbind(model=rownames(i),i) } )
names(aicVals)[1] <- "Ygen"
write.csv(aicVals, "AICvalues.csv", row.names=F)

#Get parameter estimates for best models
coefs=lapply(seq_along(mods),function(i) {
  df=dlply(metamaster.reduced,"Ygen")[[i]]
  df=df[!is.na(df$value) & !is.infinite(df$value),]
  df=groupedData(value~richness|Ref,data=df)
  bestmod.name=names(which.min(lapply(mods[[i]],AIC)))
  best.mod=mods[[i]][[bestmod.name]]
  #update(best.mod,data=df,start=c(a=fixef(best.mod)[1],b=fixef(best.mod)[2]),method="REML",verbose=T,
  #       control=nlmeControl(maxIter=1000))
  summary(best.mod)$tTable
} )
names(coefs)=names(mods); 


coefs <- ldply(coefs)
coefs$coef <- c("a", "b")
names(coefs)[1] <- "Ygen"
write.csv(coefs, "COEFvalues.csv", row.names=F)

#PLOT ALL THE THINGS
lapply(seq_along(mods), function(i){
  pdf(paste(names(mods)[[i]], ".pdf", sep=""), width=10, height=10)
  modset <- mods[[i]]
  
  #get some data for plotting
  df=dlply(metamaster.reduced,"Ygen")[[i]]
  df=df[!is.na(df$value) & !is.infinite(df$value),]
  df=groupedData(value~richness|Ref,data=df)
  
  #which mod is it
  bestmod.name <- names(which.min(lapply(modset, AIC)))
  
  #Get the mod
  mod <- modset[[bestmod.name]]
  
  #start a plot
  p <- ggplot(data=df, aes(x=richness, y=value)) + 
    geom_point(color="lightgrey", alpha=0.5) +
    theme_bw()
  lms <- quantile(df$value, c(0.005, 0.995))
  if(lms[1]>0) lms[1] <- 0
  
  p <- p+ylim(c(lms[1], lms[2]))
  
  x <- seq(0.01, max(df$richness), .01)
  if(bestmod.name=="Logarithmic"){
    y <- apply(coef(mod), 1, function(vals) vals[1]+vals[2]*log(x))  
    f <- data.frame(richness=x, value=fixef(mod)[1] + fixef(mod)[2]*log(x))
  }
  if(bestmod.name=="Power"){
    y <- apply(coef(mod), 1, function(vals) vals[1]*x^vals[2])  
    f <- data.frame(richness=x, value=fixef(mod)[1] *x^fixef(mod)[2])
  }
  
  #reorganize disaggregated data
  y <- as.data.frame(y)
  y$richness <- x
  y <- melt(y, "richness")
  
  p <- p+
    geom_line(data=y, mapping=aes(group=variable), color="lightgrey") +
    geom_line(data=f, color="black", lwd=2)
  
  
  print(p)
  
  dev.off()
})









