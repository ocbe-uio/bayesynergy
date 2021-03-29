#' Summary function for bayesynergy object
#' 
#' @description 
#' A function summarizing posterior inference of a \code{bayesynergy} object.
#' 
#' @param object The fitted model to summarise, from \code{\link{bayesynergy}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples 
#' \dontrun{
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
#' x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
#' fit <- bayesynergy(y_mat,x_mat)
#' summary(fit)
#' }
#' 
#' 
#' @export

summary.bayesynergy <- function(object,...){
  posterior = rstan::extract(object$stanfit)
  coef_names = setdiff(names(posterior),c("C","p0","p01","p02","Delta","CPO","lp__",
                                          "ec50_1", "ec50_2", "b1", "b2", 
                                          "theta_1", "theta_2", "z",
                                          "log_rVUS_f","log_rVUS_p0", "log_rVUS_Delta", "log_rVUS_syn", "log_rVUS_ant",
                                          "s2_log10_ec50_1", "s2_log10_ec50_2"))
  
  summ = rstan::summary(object$stanfit,pars=coef_names,probs=c(0.025,.5,0.975))$summary
  print(summ,digits=3)
  
  cat("\n")
  cat("log-Pseudo Marginal Likelihood (LPML) = ", object$LPML, "\n")
}
