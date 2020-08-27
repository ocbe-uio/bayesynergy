#' Summary function for bayesynergy object
#' 
#' @description 
#' A function summarizing posterior inference of a \code{bayesynergy} object.
#' 
#' @param object The fitted model to summarise, from \code{\link{bayesynergy}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples 
#' library(BayeSyneRgy)
#' data("mathews_DLBCL")
#' y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
#' x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
#' fit <- bayesynergy(y_mat,x_mat)
#' summary(fit)
#' 
#' 
#' @export

summary.bayesynergy <- function(object,...){
  posterior = rstan::extract(object$stanfit)
  coef_names = setdiff(names(posterior),c("pij_0","pij_01","pij_02","Delta_ij","CPO","lp__"))
  
  summ = rstan::summary(object$stanfit,pars=coef_names,probs=c(0.025,.5,0.975))$summary
  print(summ,digits=3)
  
  cat("\n")
  cat("log-Pseudo Marginal Likelihood (LPML) = ", object$LPML, "\n")
}
