### trying to adapt Jack's code for caterpillar plots


#### start function ####
#### Arguments:  DATASET=data frame used for building model
####             Model.out=model created by lmer function

resids.table<-function(DATASET,Model.out)
{
  #get list of species names
  DATASET<-DATASET[order(DATASET$species),]
  species.list<-unique(sort(DATASET$species))
  #construct the matrix Omega
  #tau00<-as.numeric((show(VarCorr(Model.out)))[1,3])
  #tau11<-as.numeric((show(VarCorr(Model.out)))[2,3])
  #tau22<-as.numeric((show(VarCorr(Model.out)))[3,3])
  #rho01<-as.numeric((show(VarCorr(Model.out)))[2,5])
  #rho02<-as.numeric((show(VarCorr(Model.out)))[3,5])
  #rho12<-as.numeric((show(VarCorr(Model.out)))[3,6])
  #tau01<-sqrt(tau00*tau11)*rho01
  #tau02<-sqrt(tau00*tau22)*rho02
  #tau12<-sqrt(tau22*tau11)*rho12
  tau00<-VarCorr(Expq4)$species[1,1]
  tau11<-VarCorr(Expq4)$species[2,2]
  tau22<-VarCorr(Expq4)$species[3,3]
  tau01<-VarCorr(Expq4)$species[2,1]
  tau02<-VarCorr(Expq4)$species[3,1]
  tau12<-VarCorr(Expq4)$species[3,2]
  
  
  omegaH<-matrix(c(tau00,tau01,tau02,tau01,tau11,tau12, tau02,tau12,tau22), ncol=3)
  
  #omegaH<-VarCorr(Expq4)$species
  #obtain the level-2 residuals from the model
  reslev2.u0<-coef(Model.out)$species[,1]-fixef(Model.out)[1]
  reslev2.u1<-coef(Model.out)$species[,2]-fixef(Model.out)[2]
  reslev2.u2<-coef(Model.out)$species[,3]-fixef(Model.out)[3]
  
  
  #limma is needed for the blockDiag function
  library(limma)
  
  #locate where each level-2 unit starts and stops
  cumsum(table(DATASET$species))->stops
  starts<-c(1,stops+1)
  starts<-starts[-length(starts)]
  
  #initialize the matrices V, R, and S
  Vmat<-NULL
  Rmat<-NULL
  sh<-NULL
  
  #build the block diagonal matrices species by species
  for (i in 1:length(stops)) {
    #next line creates Zh--currently an intercept and lntemp-2.8
    Zmat<-cbind(rep(1,stops[i]-starts[i]+1),
                DATASET[starts[i]:stops[i],"lntemp"]-log(15),
                (DATASET[starts[i]:stops[i],"lntemp"]-log(15))^2)  
    Rhat<-Zmat%*%omegaH
    Vpart<-Zmat%*%omegaH%*%t(Zmat)  
    Vmat<-if(i==1) Vpart else blockDiag(Vmat,Vpart)
    Rmat<-if(i==1) Rhat else blockDiag(Rmat,Rhat)
    sh<-if(i==1) omegaH else blockDiag(sh,omegaH)
  }   
  #add sigma2 to the diagonal entries
  # was: sigma2<-as.numeric((show(VarCorr(Model.out)))[4,3])
  sigma2<-(attr(VarCorr(Expq4), "sc"))^2
  V.true<-Vmat+diag(rep(sigma2,dim(Vmat)[1]))
  
  #need X matrix of model variables. Currently these are
  #intercept, lntemp-2.8, and (lntemp-2.8)^2
  xmat<-cbind(rep(1,dim(DATASET)[1]),
              DATASET[,"lntemp"]-log(15),(DATASET[,"lntemp"]-log(15))^2)
  
  #sand.wich is the middle part of variance expression—the
  #part in the parentheses plus premultiplication by V^-1
  sand.wich<-(solve(V.true))%*%(V.true-xmat%*%
                                  (solve(t(xmat)%*%(solve(V.true))%*%xmat))%*%t(xmat))
  res2.sd<-sqrt(diag(sh-t(Rmat)%*%sand.wich%*%(solve(V.true))%*%Rmat))
  
  #pull out the different random effects
  sd.u0<-res2.sd[names(res2.sd)=='1']
  sd.u1<-res2.sd[names(res2.sd)=='2']
  sd.u2<-res2.sd[names(res2.sd)=='3']
  
  
  #construct data frame of results
  out.resids<-data.frame(as.numeric(species.list),reslev2.u0,
                         reslev2.u1,reslev2.u2,sd.u0,sd.u1,sd.u2)
  colnames(out.resids)<-c("Species.num","u0", "u1", "u2", "sd.u0","sd.u1", "sd.u2")
  rownames(out.resids)<-species.list
  out.resids
}
out.resids<-resids.table(PLD,Expq4) 


cat.plot3<-function(outfile,uval) { 
  oldmar<-par('mar') 
  #expand margins 
  par(mai=c(1,1.1,.3,.3)) 
  species.list<-rownames(outfile) 
  #select correct random effects 
  u.resids<-outfile[,uval] 
  sd.uresids<-outfile[,paste("sd",uval,sep=".")] 
  #create label for y-axis 
  number<-substr(uval,2,2) 
  label<-paste(number, "i",sep="") 
  #obtain limits for 95% confidence intervals 
  u1.upper<- u.resids +qnorm(.975)*sd.uresids 
  u1.lower<- u.resids +qnorm(.025)*sd.uresids 
  plot(outfile[order(u.resids),uval], axes=FALSE, pch=17, cex=1.5, cex.lab=1.9,
       ylim=c(-4,3), xlim=c(0,75), 
       xlab='Rank order of Residuals', 
       ylab= substitute(paste('Species-level Residual ',
                              hat(u)[val]), list=list(val=label))) 
  sort.lower<-u1.lower[order(u.resids)] 
  sort.upper<-u1.upper[order(u.resids)] 
  #identify points whose confidence intervals do not overlap 0 
  #outlier.flag=1 for deviant points 
  outlier.flag<-ifelse(sort.lower>0 | sort.upper<0,1,0) 
  #use different color for error bars of deviant residuals 
  for(i in 1:(length(u.resids))) 
    if (outlier.flag[i]==0) arrows(i,sort.lower[i],i,
                                   sort.upper[i], 
                                   angle=90,code=3,length=.05,col='gray60') else arrows(i,sort.lower[i],i,sort.upper[i], 
                                                                                        angle=90,code=3,length=.05,col='tomato') 
  #abline(h=0,col=4,lty=2,lwd=2) 
  x<-cbind(0,75)
  y<-cbind(0,0)
  lines(x,y, lty=2, lwd=2, col=4)
  outfile[order(u.resids),]->sorted.outfile 
  #use different color for predicted values of deviant residuals 
  #normal points 
  points((1:length(species.list))[outlier.flag ==0], sorted.outfile[outlier.flag==0,uval],pch=17,cex=1.5) 
  #deviant points 
  points((1:length(species.list))[outlier.flag ==1], sorted.outfile[outlier.flag==1,uval],pch=17,cex=1.5,col=2) 
  #expand print area for labels 
  #par(xpd=TRUE) 
  #identify(1:length(species.list),outfile[order(u.resids),uval], 
  #rownames(outfile[order(u.resids),]),cex=.75,font=4,col=2) 
  #par(xpd=FALSE) 
  par(mar=oldmar) 
  #axis(1, c(0,25,50,75), cex.axis = 1.8, pos = min(u1.lower))
  #axis(2, c(min(u1.lower),-2,0,2,max(u1.upper)), cex.axis=1.8, las=1, pos=0)
} 

cat.plot3(out.resids, 'u1')
axis(1, c(0,25,50,75), cex.axis = 1.8, pos = -2)
axis(2, c(-2,-1,0,1,2), cex.axis=1.8, las=1, pos=0)
legend(71,1.9,'A', cex=1.7, bty='n')
