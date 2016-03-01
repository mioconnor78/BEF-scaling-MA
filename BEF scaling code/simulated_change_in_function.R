######################################################################
#' @title Simulations of change in diversity leading to changes in function
#' 
#' @author Jarrett Byrnes
#' @author jarrett.byrnes@umb.edi
#' 
#' @log
#'  2/25/2016 - First draft
######################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

#Vector of scaling coefficients
b <- c(0.25, 0.47, 0.53)


## below: consider changing sd_change, maybe the mean values of change. what are the units on those? log of % rate of change. log(S1/S2) = rate*dur.   SI = S.o* exp(RD); make sure the model structure is ok (lines 37 below)

#list of richness changes
change <- list(loss = -0.01, hold = 0, gain = 0.01)
sd_change = 0.02 

#some other params - might want to kick nsims up? Or not - it's a lot
nYears<- 20
nsims = 500

#make that simulation data frame
simDF <- data.frame(expand.grid(trophic_group = rep(1:3, nsims),
                                 scenario=c("loss", "hold", "gain")),
                     
                     s0 = round(runif(9*nsims, 5, 40))) %>%
  #calculate change in diversity
  group_by(scenario) %>%
  mutate(s1 = s0*exp(nYears*rnorm(3*nsims, change[[scenario[1]]], sd_change))) %>%
  ungroup() %>%
  
  #caculate function for each trophic group before and after
  group_by(trophic_group) %>%
  mutate(f0 = s0^b[trophic_group[1]],
         f1 = s1^b[trophic_group[1]]) %>%
  ungroup() %>%
  mutate(lr_s = log(s1)-log(s0),
         lr_f = log(f1)-log(f0))
  

#Plotting
#richness before and after
ggplot(data=simDF, aes(x=s0, y=s1))+ 
  geom_point() + 
  facet_wrap(~scenario, scale="free_y") + 
  geom_abline(slope=1, intercept=0, col="red")

#function before and after
ggplot(data=simDF, aes(x=f0, y=f1))+ 
  geom_point() +
  facet_grid(trophic_group~scenario, scale="free") + 
  geom_abline(slope=1, intercept=0, col="red")


#log ratio of richness change
ggplot(data=simDF, aes(x=lr_s))+ 
  geom_histogram() +
  facet_grid(trophic_group~scenario, scale="free") + 
  geom_vline(xintercept=0, col="red")


#log ratio of function change
ggplot(data=simDF, aes(x=lr_f))+ 
  geom_histogram() +
  facet_grid(trophic_group~scenario, scale="free") + 
  geom_vline(xintercept=0, col="red")