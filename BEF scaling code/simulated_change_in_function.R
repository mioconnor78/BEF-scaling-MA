######################################################################
#' @title Simulations of change in diversity leading to changes in function
#' 
#' @author Jarrett Byrnes
#' @author jarrett.byrnes@umb.edu
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
# note on 3/1: the loss is the rate of loss. so 1 species / decade * decade gives us log(S2/S1) per decade. So I think this is analagous to the values in marks histogram... but not the spreadsheet. Need good estimates of deltaSR/decade for this. 
# right now, the values of 0.1 are equivalent to 1 species per decade.

#list of richness changes
change <- list(loss = -0.01, hold = 0, gain = 0.01)
sd_change = 0.2

#some other params - might want to kick nsims up? Or not - it's a lot
nYears<- 20
nsims = 500

#make that simulation data frame
simDF <- data.frame(expand.grid(trophic_group = rep(1:3, nsims),
                                 scenario=c("loss", "hold", "gain")),
                     
                   # s0 = round(runif(9*nsims, 5, 40))) %>%
  s0 = round(rlnorm(9*nsims, 3, 0.83) )) %>% #from distributions of S0 in dornelas & vellend
  #calculate change in diversity
  group_by(scenario) %>%
  #  mutate(s1 = s0*exp(nYears*rnorm(3*nsims, change[[scenario[1]]], sd_change))) %>%
    mutate(s1 = s0*rnorm(3*nsims, exp(nYears*change[[scenario[1]]]), 0.2)) %>%
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
  # geom_vline(xintercept=mean(simDF$lr_s), col = 3) ## help? how do we add a line for the mean of each distribution in each panel?


#log ratio of function change
ggplot(data=simDF, aes(x=lr_f))+ 
  geom_histogram() +
  facet_grid(trophic_group~scenario, scale="free") + 
  geom_vline(xintercept=0, col="red")