## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load_example,cache=T-----------------------------------------------------
library(bayesynergy)
data("mathews_DLBCL")
y = mathews_DLBCL[[1]][[1]]
x = mathews_DLBCL[[1]][[2]]
head(cbind(y,x))

## ----regular_fit, cache = T, warning = FALSE----------------------------------
fit = bayesynergy(y,x, drug_names = c("ibrutinib", "ispinesib"),
                  units = c("nM","nM"),bayes_factor = T)

## ----summary_fit,cache=T------------------------------------------------------
summary(fit)

