library(dplyr)
dorn <- read.csv("1248484_s1.csv") %>% 
  arrange(Year) %>%
  group_by(ID) %>%
  slice(1L) %>%
  ungroup

vel <- read.csv("/Users/maryo/Documents/projects/BEF synthesis/BEF-scaling-MA/data/sd01.csv", header=T)
vel <- subset(vel, vel$SR_analysis==1) %>%
  filter(!is.na(SR_Year1_CT))


library(MASS)

par(mfrow=c(1,2))

poisVel <- fitdistr(dorn$S, "poisson")
lDorn <- fitdistr(dorn$S, "lognormal")
plot(density(dorn$S), main = "Dornelas et al.\nInitial Richness")
#matplot(0:1500, dpois(0:1500, poisDorn$estimate), add=T, col="red", type="l")
matplot(0:1500, dlnorm(0:1500, lDorn$estimate[1], lDorn$estimate[2]), add=T, col="blue", type="l")



poisVel <- fitdistr(vel$SR_Year1_CT, "poisson")
lVel <- fitdistr(vel$SR_Year1_CT, "lognormal")
plot(density(vel$SR_Year1_CT), main = "Vellend et al.\nInitial Richness")
#matplot(0:1500, dpois(0:1500, poisVel$estimate), add=T, col="red", type="l")
matplot(0:1500, dlnorm(0:1500, lVel$estimate[1], lDorn$estimate[2]), add=T, col="blue", type="l")

par(mfrow=c(1,1))
