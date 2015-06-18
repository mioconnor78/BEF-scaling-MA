#For Fig1
library(ggplot2)
library(dplyr)

data <- SST5 #<- read.csv("sst5_20150601_jeb.csv")

#I want some new models!
modBtrophic<-lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1 + logSc*factor(TG1) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)
modBtrophic2<-lmer(logY.rs ~ logSc+log(Tscale) + Sys1 + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

modBtrophic3 <- lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1*factor(TG1) - logSc:Sys1T:factor(TG1)2 + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

## MO use this model:
modBtrophic3 <- lmer(logY.rs ~ logSc*log(Tscale) + logSc*Sys1*factor(TG1) - logSc:Sys1:factor(TG1) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = FALSE, na.action=na.omit)

SST5$modBtrophic.fitted <- fitted(modBtrophic)
SST5$modBtrophic3.fitted <- fitted(modBtrophic3)



#####
# Make curves separated by Sys1
#####

modBtrophicPred.varyNone <- data.frame(logSc=rep(log(seq(1,5,length.out=100))-log(8), 2), Sys1=factor(c(rep("A",100), rep("T",100))))

modBtrophicPred.varyNone <- data.frame(logSc=rep(log(seq(1,5,length.out=100))-log(8), 2), TG1=factor(c(rep("1",100), rep("2",100), rep("4",100))))


modBtrophicPred.varyNone$logY.rs.predicted <-  fixef(modBtrophic3)[1] + fixef(modBtrophic3)[2]*modBtrophicPred.varyNone$logSc + mean(SST5$Tscale)*fixef(modBtrophic3)[3]

modBtrophicPred.varyNone$logY.rs.predicted[101:200] <- modBtrophicPred.varyNone$logY.rs.predicted[101:200]+fixef(modBtrophic)[4]

#FIG 1!!!!
fitted_plot + facet_grid(. ~ Sys1) +
  geom_line(data=modBtrophicPred.varyNone, mapping=aes(x=exp(logSc), y=exp(logY.rs.predicted), group=Sys1), lwd=2)

####
#Try and get TG1 in there?
####

###The master plot
fitted_plot <- ggplot() +
  geom_line(data=SST5, aes(x=exp(logSc), y=exp(modBtrophic3.fitted), group=Entry), color="lightgrey") +
  theme_bw(base_size=17)
fitted_plot

#create some predictions - the intercepts all come out too low, and I'm not sure why...
modBtrophicPred.varyall <- data.frame(expand.grid(logSc=log(seq(1,5,.01))-log(8), Tscale=mean(SST5$Tscale),
                                          Sys1=factor(levels(SST5$Sys1)), 
                                          TG1=factor(unique(SST5$TG1)))) %>%
  filter( !(Sys1=="T" & TG1==2)) %>%
  mutate(logY.rs.predicted=predict(modBtrophic, 
                                   newdata=data.frame(Sys1=Sys1, TG1=TG1, Tscale=Tscale, logSc=logSc), 
                                   re.form=~0, type="response"))


#not sure why the intercepts are not high enough here...
fitted_plot + facet_grid(Sys1 ~ TG1) +
  geom_line(data=modBtrophicPred.varyall, mapping=aes(x=exp(logSc), y=exp(logY.rs.predicted)))

#######
#And this gets a "non-conformable aguments' error. ARGH. predic.mermod, why do you hate me so
#######

modBtrophic3Pred.varyall <- data.frame(expand.grid(logSc=log(seq(1,5,.01))-log(8), Tscale=mean(SST5$Tscale),
                                                  Sys1=factor(levels(SST5$Sys1)), 
                                                  TG1=factor(unique(SST5$TG1)))) %>%
                           filter( !(Sys1=="T" & TG1=="2")) 


#Had to roll my own
makeModB3Pred <- function(newdata, object=modBtrophic3){
  X_orig <- getME(object, "X")
  form <- formula(object, 
                  fixed.only = TRUE)
  
  form <- form[[length(form)]]
  
  #RHS <- formula(substitute(~R,
  #                          list(R=RHSForm(formula(object,fixed.only=TRUE)))))
  
  RHS <- formula(substitute(~R,
                            list(R=form)))
  
  
  Terms <- terms(object, fixed.only = TRUE)
  isFac <- vapply(mf <- model.frame(object, fixed.only = TRUE), 
                  is, "factor", FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) NULL else lapply(mf[isFac], levels)
  
  mfnew <- model.frame(delete.response(Terms), 
                       newdata, na.action = na.pass, xlev = orig_levs)
  
  X <- model.matrix(RHS, mfnew, contrasts.arg = attr(X_orig, 
                                                     "contrasts"))
  X <- X[,-c(11,13)]
  
  offset <- rep(0, nrow(X))
  tt <- terms(object)
  
  fit.na.action <- attr(mfnew, "na.action")
  
  
  pred <- drop(X %*% fixef(object))
  pred <- pred + offset
  
  pred
}


modBtrophic3Pred.varyall$logY.rs.predicted <- makeModB3Pred(modBtrophic3Pred.varyall)

mean(rnorm(50))

50 %>% 
  rnorm %>% 
  mean

fitted_plot3 <- SST5 %>% 
  mutate(Sys1 = as.character(Sys1),
         Sys1 = ifelse(Sys1 == "T", "Terrestrial", "Aquatic")) %>%
  ggplot(aes(x=exp(logSc),
             y=exp(modBtrophic3.fitted),
             group=Entry)) +
  geom_line(color="black", alpha = "0.1", size = 1.5) +
  theme_bw(base_size=17) +
  facet_grid(Sys1 ~ TG1, scales = "free_y") +
  ylab("trophic level") +
  xlab("species richness")

fitted_plot3+
  geom_line(data=modBtrophic3Pred.varyall, mapping=aes(x=exp(logSc), y=exp(logY.rs.predicted), group = NULL), colour = "red", size = 2)

ggsave("bdef.png")
getwd()
