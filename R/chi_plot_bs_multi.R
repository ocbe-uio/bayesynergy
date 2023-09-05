# synergy screen 

#' data("mathews_DLBCL")
experiment1 <- list(y = mathews_DLBCL[[1]][[1]],
                    x = mathews_DLBCL[[1]][[2]],
                    drug_names = c("ispinesib","ibrutinib"))
experiment2 <- list(y = mathews_DLBCL[[2]][[1]],
                    x = mathews_DLBCL[[2]][[2]],
                    drug_names = c("canertinib","ibrutinib"))
experiments <- list(experiment1,experiment2)

# fit <- synergyscreen(experiments)
# plot(fit)

x <- fit


plot_synergy_antagonism <- function(fitobj){
  
  
  
  
}



process_ss <- function(fitobj){
  
  
  
  
}



