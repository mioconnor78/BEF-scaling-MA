# simple plot
library(sjPlot)
library(arm)

#0# choose the model to use
data <- SST5

mod4 <- lmer(logY.rs ~ logSc*Sys1*TG1 + log(Tscale) + (1 + logSc|Entry) + (1 + logSc|ExptA) + (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod <- mod4

mod2F <- lmer(logY.rs ~ logSc + log(Tscale) + (1 + logSc|Entry) +  (1 + logSc|ExptA) +  (1 + logSc|Study), data=data, REML = TRUE, na.action=na.omit)

mod <- mod2F


pdf(file = "catplotEntry.pdf", width = 7.5, height = 3)
# make the caterpillar plot for Entry level coefs
cat.plot <- sjp.lmer(mod, sort = "logSc", ri.nr = 1, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

pdf(file = "catplotExptA.pdf", width = 7.5, height = 3)
# make the caterpillar plot for ExptA level coefs
cat.plot2 <- sjp.lmer(mod, sort = "logSc", ri.nr = 2, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

pdf(file = "catplotStudy.pdf", width = 7.5, height = 3)
# make the caterpillar plot for Study level coefs
cat.plot3 <- sjp.lmer(mod, sort = "logSc", ri.nr = 3, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

