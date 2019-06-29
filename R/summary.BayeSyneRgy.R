#' Summary function for BayeSyneRgy object
#' 
#' @description 
#' A function summarizing posterior inference of a \code{BayeSyneRgy} object.
#' 
#' @param object The fitted model to summarise, from \code{\link{BayeSyneRgy}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' 
#' 
#' @export

summary.BayeSyneRgy <- function(object,...){
  
  Summary_Output <- object$Summary_Output
  
  LPML <- object$OUTPUT$MCMC_Output$LPML
  S2_EPS <- object$OUTPUT$MCMC_Output$S2_EPS
  
  n_save <- length(Summary_Output[[1]])
  summary_out <- (cbind(sapply(Summary_Output, mean),
                        sapply(Summary_Output, sd),
                        sapply(Summary_Output, quantile, probs = 0.025),
                        sapply(Summary_Output, quantile, probs = 0.5),
                        sapply(Summary_Output, quantile, probs = 0.975),
                        sapply(Summary_Output, function(x,n){effectiveSize(as.mcmc(x))/n}, n = n_save)))
  colnames(summary_out) <- c("mean", "sd", paste(2.5, "%"), paste(50, "%"), paste(97.5, "%"), "ESS/N")
  print(summary_out,digits=2)
  
  cat("\n")
  cat("log-Pseudo Marginal Likelihood (LPML) = ", LPML, "\n")
  
  #Check Bio-metrics
  cat("\n")
  if(sd(S2_EPS) > 1){
    cat("S2_EPS has std. dev. > 1", "\n")
  }
  
}
