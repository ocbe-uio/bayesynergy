## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=F,cache=T-----------------------------------------------------------
df = data.frame("BF10" = c("1 to 3.2","3.2 to 10","10 to 100",">100"),"Evidence against M0" = c("Not worth more than a bare mention","Substantial","Strong","Decisive"))
knitr::kable(df,col.names = c("$\\text{BF}_{10}$","Evidence against $\\mathcal{M}_0$"),align="cc",format = "simple")

