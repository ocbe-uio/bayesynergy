#' Summary function for bayesynergy3 object
#' 
#' @description 
#' A function summarizing posterior inference of a \code{bayesynergy3} object.
#' 
#' @param object The fitted model to summarise, from \code{\link{bayesynergy3}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' 
#' @export

summary.bayesynergy3 <- function(object,...){
  posterior = rstan::extract(object$stanfit)
  coef_names = setdiff(names(posterior),c("C","p0_12","p0_13","p0_23","p0_123",
                                          "p01","p02","p03",
                                          "Delta_12","Delta_13","Delta_23","Delta_123",
                                          "CPO","GP","lp__",
                                          "ec50_1", "ec50_2", "ec50_3", "b1", "b2", 
                                          "theta_1", "theta_2","theta_3", 
                                          "z_12", "z_13", "z_23", "z_123",
                                          "log_rVUS_f","log_rVUS_p0", "log_rVUS_Delta", "log_rVUS_syn", "log_rVUS_ant",
                                          "s2_log10_ec50_1", "s2_log10_ec50_2", "s2_log10_ec50_3",
                                          "log10_ec50_1_raw","log10_ec50_2_raw"))
  
  summ = rstan::summary(object$stanfit,pars=coef_names,probs=c(0.025,.5,0.975))$summary
  print(summ,digits=3)
  
  # cat("\n")
  # cat("log-Pseudo Marginal Likelihood (LPML) = ", object$LPML, "\n")
  # if (object$model$bayes_factor){
    # cat("Estimated Bayes factor in favor of full model over non-interaction only model: ", object$bayesfactor, "\n")
  # }
}
