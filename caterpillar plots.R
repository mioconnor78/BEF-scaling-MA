# simple plot
library(sjPlot)
library(arm)

mod <- modBtrophic
mod <- modBasic
mod <- modBtrophic1

sjp.lmer(mod, sort = "logSc", fade.ns = TRUE, free.scale = TRUE, geom.colors = c(1, 1), showValueLabels = FALSE)
sjp.lmer(mod, type = "fe.cor")
sjp.lmer(mod, type = "re.qq")
