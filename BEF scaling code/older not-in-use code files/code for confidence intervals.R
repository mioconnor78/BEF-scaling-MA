### trying this script from the web to make confidence intervals on lmer outputs. 

### http://rpubs.com/hughes/116374

library(fda)

#####LM example######
#we measured plant biomass for 120 pots under 3 nutrient treatments and across a gradient of CO2
#due to limited place in our greenhouse chambers we had to use 4 of them, so we established a blocking design
data<-data.frame(C=runif(120,-2,2),N=gl(n = 3,k = 40,labels = c("Few","Medium","A_lot")),Block=rep(rep(paste0("B",1:4),each=10),times=3)) 

xtabs(~N+Block,data)

modmat <- model.matrix(~Block+C*N,data)  #this is cool; I don't know why he did this

#the paramters of the models
params<-c(10,-0.4,2.3,-1.5,1,0.5,2.3,0.6,2.7)

#simulate a response vector
data$Biom<-rnorm(120,modmat%*%params,1)

#fit the model
m<-lm(Biom~Block+C*N,data)
summary(m)

## so far, this is vary analagous to what I try to do, with different levels of N that I want to plot along C.

# Here we would normally continue and make some model checks. As output from the model we would like to plot the effect of CO2 on plant biomass for each level of N addition. Of course we want to average out the Block effect (otherwise we would have to plot one separate line for each Block). This is how it works:
 
#his new function:  
plot_fit<-function(m,focal_var,inter_var=NULL,RE=NULL,n=20,n_core=4){
  require(arm)  
  dat<-model.frame(m)
  #turn all character variable to factor
  dat<-as.data.frame(lapply(dat,function(x){
    if(is.character(x)){
      as.factor(x)
    }
    else{x}
  }))
  #make a sequence from the focal variable
  x1<-list(seq(min(dat[,focal_var]),max(dat[,focal_var]),length=n))
  #grab the names and unique values of the interacting variables
  isInter<-which(names(dat)%in%inter_var)
  if(length(isInter)==1){
    x2<-list(unique(dat[,isInter]))
    names(x2)<-inter_var
  }
  if(length(isInter)>1){
    x2<-lapply(dat[,isInter],unique)
  }
  if(length(isInter)==0){
    x2<-NULL
  }
  #all_var<-x1
  #add the focal variable to this list
  all_var<-c(x1,x2)
  #expand.grid on it
  names(all_var)[1]<-focal_var
  all_var<-expand.grid(all_var)
  
  #remove varying variables and non-predictors
  dat_red<-dat[,-c(1,which(names(dat)%in%c(focal_var,inter_var,RE,"X.weights."))),drop=FALSE]
  if(dim(dat_red)[2]==0){
    new_dat<-all_var
  }
  else{
    fixed<-lapply(dat_red,function(x) if(is.numeric(x)) mean(x) else factor(levels(x)[1],levels = levels(x)))
    #the number of rows in the new_dat frame
    fixed<-lapply(fixed,rep,dim(all_var)[1])
    #create the new_dat frame starting with the varying focal variable and potential interactions
    new_dat<-cbind(all_var,as.data.frame(fixed)) 
    #get the name of the variable to average over, debug for conditions where no variables are to be avergaed over
    name_f<-names(dat_red)[sapply(dat_red,function(x) ifelse(is.factor(x),TRUE,FALSE))]
  }  
  
  
  #get the predicted values
  cl<-class(m)[1]
  if(cl=="lm"){
    pred<-predict(m,newdata = new_dat,se.fit=TRUE)
  }
  
  if(cl=="glm" | cl=="negbin"){
    #predicted values on the link scale
    pred<-predict(m,newdata=new_dat,type="link",se.fit=TRUE)
  }
  if(cl=="glmerMod" | cl=="lmerMod"){
    pred<-list(fit=predict(m,newdata=new_dat,type="link",re.form=~0))
    #for bootstrapped CI
    new_dat<-cbind(new_dat,rep(0,dim(new_dat)[1]))
    names(new_dat)[dim(new_dat)[2]]<-as.character(formula(m)[[2]])
    mm<-model.matrix(formula(m,fixed.only=TRUE),new_dat)
  }
  #average over potential categorical variables  
  if(length(name_f)>0){
    if(cl=="glmerMod" | cl=="lmerMod"){
      coef_f<-lapply(name_f,function(x) fixef(m)[grep(paste0("^",x,"\\w+$"),names(fixef(m)))])
    }
    else{
      coef_f<-lapply(name_f,function(x) coef(m)[grep(paste0("^",x,"\\w+$"),names(coef(m)))])
    }    
    pred$fit<-pred$fit+sum(unlist(lapply(coef_f,function(x) mean(c(0,x)))))
  }
  #to get the back-transform values get the inverse link function
  linkinv<-family(m)$linkinv
  
  #get the back transformed prediction together with the 95% CI for LM and GLM
  if(cl=="glm" | cl=="lm"){
    pred$pred<-linkinv(pred$fit)
    pred$LC<-linkinv(pred$fit-1.96*pred$se.fit)
    pred$UC<-linkinv(pred$fit+1.96*pred$se.fit)
  }
  
  #for GLMM need to use bootstrapped CI, see ?predict.merMod
  if(cl=="glmerMod" | cl=="lmerMod"){
    pred$pred<-linkinv(pred$fit)
    predFun<-function(.) mm%*%fixef(.)
    bb<-bootMer(m,FUN=predFun,nsim=200,parallel="multicore",ncpus=n_core) #do this 200 times
    bb$t<-apply(bb$t,1,function(x) linkinv(x))
    #as we did this 200 times the 95% CI will be bordered by the 5th and 195th value
    bb_se<-apply(bb$t,1,function(x) x[order(x)][c(5,195)])
    pred$LC<-bb_se[1,]
    pred$UC<-bb_se[2,] 
  }
  
  #the output
  out<-as.data.frame(cbind(new_dat[,1:(length(inter_var)+1)],pred$LC,pred$pred,pred$UC))
  names(out)<-c(names(new_dat)[1:(length(inter_var)+1)],"LC","Pred","UC")
  return(out)
}




pred <- plot_fit(m, focal_var = "C", inter_var = "N")
head(pred)




###### switching approach to the recommendation at the bottom of http://glmm.wikidot.com/faq

### estimate Confidence intervals on conditional means/BLUPs/random effects

library("lme4")
library("ggplot2") # Plotting

data("Orthodont", package="MEMSS")

fm1 <- lmer(
  formula = distance ~ age*Sex + (age|Subject)
  , data = Orthodont
)

newdat <- expand.grid(
  age=c(8,10,12,14)
  , Sex=c("Female","Male")
  , distance = 0
)

mm <- model.matrix(terms(fm1),newdat)

newdat$distance <- predict(fm1,newdat,re.form=NA)

## or newdat$distance <- mm %*% fixef(fm1)

pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))
tvar1 <- pvar1+VarCorr(fm1)$Subject[1]  ## must be adapted for more complex models
cmult <- 2 ## could use 1.96
newdat <- data.frame(
  newdat
  , plo = newdat$distance-cmult*sqrt(pvar1)
  , phi = newdat$distance+cmult*sqrt(pvar1)
  , tlo = newdat$distance-cmult*sqrt(tvar1)
  , thi = newdat$distance+cmult*sqrt(tvar1)
)

#plot confidence
library(ggplot2)
g0 <- ggplot(newdat, aes(x=age, y=distance, colour=Sex))+geom_point()
g0 + geom_errorbar(aes(ymin = plo, ymax = phi))+
  labs(title="CI based on fixed-effects uncertainty ONLY")
#plot prediction
g0 + geom_errorbar(aes(ymin = tlo, ymax = thi))+
  labs(title="CI based on FE uncertainty + RE variance")


##### 
## building on the above, from http://rpubs.com/hughes/116374
#### 

#first case simple lmer, simulate 100 data points from 10 groups with one continuous fixed effect variable
x<-runif(100,0,10)
f<-gl(n = 10,k = 10)
data<-data.frame(x=x,f=f)
modmat<-model.matrix(~x,data)
#the fixed effect coefficient
fixed<-c(1,0.5)
#the random effect
rnd<-rnorm(10,0,0.7)
#the simulated response values
data$y<-rnorm(100,modmat%*%fixed+rnd[f],0.3)

#model
m<-lmer(y~x+(1|f),data)

#second version with bootMer
#we have to define a function that will be applied to the nsim simulations
#here we basically get a merMod object and return the fitted values
predFun<-function(.) mm %*% fixef(.) 
bb<-bootMer(m,FUN=predFun,nsim=200) #do this 200 times
#as we did this 200 times the 95% CI will be bordered by the 5th and 195th value
bb_se<-apply(bb$t,2,function(x) x[order(x)][c(5,195)])
newdat$blo<-bb_se[1,]
newdat$bhi<-bb_se[2,]

plot(y~x,data)
lines(newdat$x,newdat$y,col="red",lty=2,lwd=3)
lines(newdat$x,newdat$plo,col="blue",lty=2,lwd=2)
lines(newdat$x,newdat$phi,col="blue",lty=2,lwd=2)
lines(newdat$x,newdat$tlo,col="orange",lty=2,lwd=2)
lines(newdat$x,newdat$thi,col="orange",lty=2,lwd=2)
lines(newdat$x,newdat$bhi,col="darkgreen",lty=2,lwd=2)
lines(newdat$x,newdat$blo,col="darkgreen",lty=2,lwd=2)
legend("topleft",legend=c("Fitted line","Confidence interval","Prediction interval","Bootstrapped CI"),col=c("red","blue","orange","darkgreen"),lty=2,lwd=2,bty="n")

