#' Bayesian semi-parametric modelling for in-vitro drug combination experiments
#'
#' @description
#' The function \code{bayesynergy3} is the three drug extension of \code{\link[bayesynergy:bayesynergy]{bayesynergy}}. It will fit a Bayesian semi-parametric model for in-vitro drug combination experiments to estimate synergistic and antagonistic effects for three drug combinations.
#'
#' @param y vector or matrix of viability measures. Replicates can be given in long or wide format.
#' @param x three-column matrix of drug concentrations.
#' @param type integer; the type of model used. Must be one of the following: 1 (Splines), 2 (GP with squared exponential kernel), 3 (GP with Matérn kernel) or 4 (GP with rational quadratic kernel).
#' @param drug_names vector of size 2; names of the drugs utilized for the experiment.
#' @param experiment_ID character; identifier of experiment, typically name of cell Line.
#' @param units vector of size 2; concentration units for the drugs, e.g. c("\eqn{\mu}M","\eqn{\mu}M")
#' @param lower_asymptotes logical; if TRUE the model will estimate the lower asymptotes of monotherapy curves.
#' @param lambda numeric; the parameter controls the residual noise observed in the heteroscedastic model when f = 0.
#' @param control list; passed on to the stan sampler, e.g. for setting adapt_delta.
#' @param ... Arguments passed to \code{\link[rstan:sampling]{rstan::sampling}} or \code{\link[rstan:vb]{rstan::vb}} (e.g. iter, chains).
#'
#' @return An object of S3 class "\code{bayesynergy3}", which is a list with the following entries
#' \tabular{ll}{
#' stanfit \tab An object of class \code{\link[rstan:stanmodel-class]{stanmodel}}, returned from the sampler. \cr
#' posterior_mean \tab A list containing the posterior means of model parameters. \cr
#' data \tab A list containing the original data used to fit the model. \cr
#' model \tab A list containing model specification for the model fit. \cr
#' }
#'
#'
#' 
#' 
#' 
#' @importFrom utils modifyList
#' @export

bayesynergy3 <- function(y, x, type = 3, drug_names=NULL, experiment_ID = NULL, units = NULL,
                         lower_asymptotes = T, lambda = 0.005, method = "sampling",
                         control = list(), ...){
  # Keep original data
  y.original = y
  x.original = x
  
  
  # Checking that data input is valid
  if (!is.numeric(x)){
    if (!(ncol(x)==3)){
      stop("Argument 'x' should be a two-column matrix of drug concentrations")
    }
  }
  if (nrow(as.matrix(y))!=nrow(x)){
    stop("Dimension mismatch! Arguments 'y' and 'x' should be equally sized.")
  }
  if (ncol(as.matrix(x))!=3){
    stop("Dimension mismatch! Argument 'x' should be a matrix with three columns containing drug concentrations.")
  }
  if (min(x[,1]) > 0 | min(x[,2]) > 0 | min(x[,3]) > 0){
    stop("Missing monotherapy data! Make sure 'x' contains zero concentrations")
  }
  # GP models need adapt_delta > 0.9, so set that here
  control.default = list(adapt_delta = 0.8)
  if (type != 1){
    control.default = list(adapt_delta = 0.9)
  }
  control = modifyList(control.default,control)
  # Setting default names for drugs and experimentID
  if (is.null(drug_names)){
    drug_names = c("Drug1", "Drug2","Drug3")
  }
  if (is.null(experiment_ID)){
    experiment_ID = "Experiment"
  }
  if (is.null(units)){
    units = c("conc.","conc.","conc.")
  }
  
  # Setting up data for the sampler
  unqX1 = log10(sort(unique(x[,1])))[-1] # Removing -Inf here
  unqX2 = log10(sort(unique(x[,2])))[-1] # Removing -Inf here
  unqX3 = log10(sort(unique(x[,3])))[-1] # Removing -Inf here
  
  n1 = length(unqX1)
  n2 = length(unqX2)
  n3 = length(unqX3)
  
  # Need to find coordinates for the observed variables in this new coordinate system
  X1 = c(0,10^unqX1)
  X2 = c(0,10^unqX2)
  X3 = c(0,10^unqX3)
  Xgrid = expand.grid(X1,X2,X3)
  Xgrid = Xgrid[order(Xgrid[,"Var3"],Xgrid[,"Var1"],Xgrid[,"Var2"]),]
  ii_obs = match(data.frame(t(round(x,digits=15))),data.frame(t(round(Xgrid,digits=15))))
  
  
  # If replicates are given in matrix form, needs to be handled by replicating coordinates.
  y = as.matrix(y)
  nrep = ncol(y)
  if (nrep > 1){
    ii_obs = rep(ii_obs,nrep) # Replicating
    ii_obs = ii_obs[which(!is.na(y))] # Removing those that are missing
    nmissing = length(rep(1:nrow(Xgrid),nrep))-length(ii_obs) # How many missing?
  } else {
    nrep = max(table(ii_obs))
    ii_obs = ii_obs[which(!is.na(y))] # Removing those that are missing
    nmissing = (n1+n2+  n3 +
                  n1*n2 + n1*n3 +
                  n2*n3 +
                  n1*n2*n3+1)*nrep-length(ii_obs)
  }
  
  y = as.vector(y)
  y = y[which(!is.na(y))] # Removing missing
  
  
  # Setting up data for Stan
  stan_data = list(n1=n1, n2=n2, n3=n3,
                   x1=unqX1, x2=unqX2, x3=unqX3, nrep=nrep,
                   y=y, nmissing=nmissing, ii_obs = ii_obs, est_la = lower_asymptotes,
                   lambda = lambda)
  # stan_data_nointeraction = stan_data
  
  # Placeholder for fit
  fit = c()
  # Placeholder for model fit without interaction term
  fitH0 = c()
  # Return code will capture warnings
  returnCode = 0
  # messages will store them
  messages = c()
  withCallingHandlers(
    warning = function(cnd){
      returnCode <<- 1
      messages <<- c(messages,cnd$message)
    },
    if (method=="sampling"){
      
      fit = rstan::sampling(stanmodels$gp_grid3,stan_data, control = control, ...)
      # if (bayes_factor){
      #   fitH0 = rstan::sampling(stanmodels$nointeraction,stan_data_nointeraction, control = control, ...)
      # }
      
    } else if (method == "vb"){
      
      fit = rstan::vb(stanmodels$gp_grid,stan_data, algorithm = "fullrank", ...)
      
    }
  )
  
  
  # # Some calculations:
  posterior = rstan::extract(fit)
  n.save = length(posterior$lp__)
  # LPML = -sum(log(apply(posterior$CPO,2,sum)/n.save))
  coef_names = names(posterior)
  # # Remove those we don't care about
  coef_names = setdiff(coef_names,c("z_12","z_13","z_23","z_123","p01","p02","p03","p0_12","p0_13","p0_23","p0_123",
                                    "Delta_12","Delta_13","Delta_23","Delta_123","wout","CPO","lp__"))
  
  # Surfaces
  # Mean response
  # First we create the arrays we need and find posterior means
  f = array(NA,dim=c(n.save,(n3+1),(n2+1),(n1+1)))
  p0 = array(NA,dim=c(n.save,(n3+1),(n2+1),(n1+1)))
  Delta = array(NA,dim=c(n.save,(n3+1),(n2+1),(n1+1)))
  # Fill in data
  # First what happens at (0,1)
  p0[,1,1,1] = 1 # At (0,1) everything is 1
  Delta[,1,1,1] = 0
  # Then monotherapies
  p0[,1,1,2:(n1+1)] = posterior$p01
  p0[,1,2:(n2+1),1] = posterior$p02
  p0[,2:(n3+1),1,1] = posterior$p03
  Delta[,1,1,2:(n1+1)] = 0
  Delta[,1,2:(n2+1),1] = 0
  Delta[,2:(n3+1),1,1] = 0
  
  # Then the two-drug effects
  p0[,1,2:(n2+1),2:(n1+1)] = posterior$p0_12
  p0[,2:(n3+1),1,2:(n1+1)] = posterior$p0_13
  p0[,2:(n3+1),2:(n2+1),1] = posterior$p0_23
  Delta[,1,2:(n2+1),2:(n1+1)] = posterior$Delta_12
  Delta[,2:(n3+1),1,2:(n1+1)] = posterior$Delta_13
  Delta[,2:(n3+1),2:(n2+1),1] = posterior$Delta_23
  # Finally the three drug effects
  p0[,2:(n3+1),2:(n2+1),2:(n1+1)] = posterior$p0_123
  Delta[,2:(n3+1),2:(n2+1),2:(n1+1)] = posterior$Delta_123
  # Construct the dose-response by summing
  f = p0 + Delta
  
  # Average over samples
  f_mean = apply(f,c(2,3,4),mean)
  p0_mean = apply(p0,c(2,3,4),mean)
  Delta_mean = apply(Delta,c(2,3,4),mean)
  
  
  
  # colnames(f_mean) = signif(c(0,10^unqX1),4)
  # rownames(f_mean) = signif(c(0,10^unqX2),4)
  
  
  # Pull out posterior mean and add some stuff
  posterior_mean = as.list(rstan::summary(fit,pars=coef_names)$summary[,'mean'])
  posterior_mean$f = f_mean
  posterior_mean$p0 = p0_mean
  posterior_mean$Delta = Delta_mean
  #
  data = list(y = y.original, x = x.original, drug_names = drug_names, experiment_ID = experiment_ID, units = units, indices = ii_obs)
  model = list(lower_asymptotes = lower_asymptotes, lambda = lambda)
  
  # Inspect residuals
  # residuals = y - as.vector(f_mean)[ii_obs]
  # # Pull out pointwise observation noise
  # point_sd = sqrt(posterior_mean$s^2*(f_mean[ii_obs]+lambda))
  # # If any residuals are outside 3*sd, user is alerted
  # if (sum(abs(residuals) > 3*point_sd) > 0){
  # res_msg = paste("Largest residuals from posterior median is",signif(max(residuals),4),", which is more than three times the observation noise. This could indicate the presence of an outlier.")
  # warning(res_msg,call. = F)
  # returnCode = 1
  # messages = c(messages,res_msg)
  # }
  #
  #
  #
  #
  object = list(stanfit = fit, posterior_mean = posterior_mean,
                data = data, model = model, returnCode = returnCode)#,
  # LPML = LPML)
  # # Set some diagnostics
  # object$divergent = NA
  # if (returnCode){
  #   object$messages = messages
  #   if (method=="sampling"){
  #     object$divergent = sum(rstan::get_divergent_iterations(fit))
  #   }
  # }
  #
  # # Calculate the Bayes Factor
  # if (bayes_factor){
  #   message("Calculating the Bayes factor")
  #   full = bridge_sampler(fit,silent=T)
  #   no_interaction = bridge_sampler(fitH0,silent=T)
  #   bfactor = bf(full,no_interaction)
  #   object$bayesfactor = bfactor$bf
  # }
  #
  #
  class(object) <- "bayesynergy3"
  return(object)
}
