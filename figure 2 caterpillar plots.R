# simple plot
library(sjPlot)
library(arm)

#0# choose the model to use
mod <- modBtrophic
mod <- modBasic


pdf(file = "catplotBtrophicEntry.pdf", width = 7.5, height = 3)
# make the caterpillar plot for Entry level coefs
cat.plot <- sjp.lmer(mod, sort = "logSc", ri.nr = 1, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

pdf(file = "catplotBtrophicExptA.pdf", width = 7.5, height = 3)
# make the caterpillar plot for ExptA level coefs
cat.plot2 <- sjp.lmer(mod, sort = "logSc", ri.nr = 2, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

pdf(file = "catplotBtrophicStudy.pdf", width = 7.5, height = 3)
# make the caterpillar plot for Study level coefs
cat.plot3 <- sjp.lmer(mod, sort = "logSc", ri.nr = 3, fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
dev.off()

#this appears to maybe be just plotting ranefs based on expt-level random effects. 

## other plots
sjp.lmer(mod, type = "fe.cor")
sjp.lmer(mod, type = "re.qq")
