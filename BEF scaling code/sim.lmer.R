library(arm)
library(ggplot2)

#random intercepts
group_int <- rnorm(10,0,100)

#create a data frame
adf <- data.frame(expand.grid(x=1:100, group = 1:10))

#y with a random intercept
adf$y <- rnorm(nrow(adf), adf$x+group_int[adf$group], 40)

#show the data visualizing the random intercept
qplot(x,y,data=adf, color=factor(group))

#model it
E1 <- lmer (y ~ x + (1 | group), data=adf)

#what did we get?
display(E1)

#simulated data to generate predictions
E1.sim <- sim(E1)

new_df <- data.frame(x=rnorm(100, 50,20), group=round(runif(100, 1,10)))

#now, different prediction types with different things included...
new_df$y_fixed <- fixef(E1.sim)[,1] +  fixef(E1.sim)[,2]*new_df$x
new_df$y_fixed_error <- rnorm(100, new_df$y_fixed, sigma.hat(E1.sim)) #could also have just gotten it from model
#new_df$y_with_ran <- new_df$y_fixed + ranef(E1)$group[new_df$group,1] #could also have done this with ranef(E1.sim) but it's complicated
new_df$y_with_ran <- new_df$y_fixed + sapply(1:100, function(i) a[[1]][i, new_df$group[1], 1]) #could also have done this with ranef(E1.sim) but it's complicated
new_df$y_with_ran_error <- rnorm(100, new_df$y_with_ran, sigma.hat(E1.sim)) #could also have done this with ranef(E1.sim) but it's complicated




