#' Bayesian semi-parametric modelling for in-vitro drug combination experiments
#' 
#' @description 
#' The function \code{bayesynergy} is the main function of the BayeSyneRgy package. It will fit a Bayesian semi-parametric model for in-vitro drug combination experiments to estimate synergistic and antagonistic effects. 
#' 


bayesynergy <- function(y_mat, x_mat, drug_names=NULL, experiment_ID = NULL, log10_conc = FALSE, lower_asymptote = T, type = 2, ...){
  object = c()
  class(object) <- "bayesynergy"
  return(object)
}
  