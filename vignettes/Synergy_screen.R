## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----synergyscreen_load, cache = T--------------------------------------------
library(bayesynergy)
data("ONeil_A375")
length(ONeil_A375$failed)

## ----synergyscreen_load2, cache = T-------------------------------------------
failed_experiment = ONeil_A375$failed[[1]]
names(failed_experiment)

## ---- echo = F,cache=T--------------------------------------------------------
colnames(failed_experiment$x) = failed_experiment$drug_names
colnames(failed_experiment$y) = "viability"

## ----cache=T------------------------------------------------------------------
head(cbind(failed_experiment$y,failed_experiment$x))

## ----synergyscreen_fit, cache = T, warning = F--------------------------------
fit_screen = synergyscreen(ONeil_A375, save_raw = F, save_plots = F, parallel = F, 
                           bayesynergy_params = list(method = "vb"))

## ----synergyscreen, cache = T, warning = FALSE, message = FALSE, fig.dim = c(8,8), fig.show="hold", results="hide", fig.keep ="all"----
plot(fit_screen)

