#' Bayesian semi-parametric modelling for in-vitro drug combination experiments
#'
#' @description
#' The function \code{bayesynergy} is the main function of the bayesynergy package. It will fit a Bayesian semi-parametric model for in-vitro drug combination experiments to estimate synergistic and antagonistic effects.
#'
#' @param y vector or matrix of viability measures. Replicates can be given in long or wide format.
#' @param x two-column matrix of drug concentrations.
#' @param type integer; the type of model used. Must be one of the following: 1 (Splines), 2 (GP with squared exponential kernel), 3 (GP with Matérn kernel) or 4 (GP with rational quadratic kernel).
#' @param drug_names vector of size 2; names of the drugs utilized for the experiment.
#' @param experiment_ID character; identifier of experiment, typically name of cell Line.
#' @param units vector of size 2; concentration units for the drugs, e.g. c("\eqn{\mu}M","\eqn{\mu}M")
#' @param lower_asymptotes logical; if TRUE the model will estimate the lower asymptotes of monotherapy curves.
#' @param heteroscedastic logical; if TRUE, the model will assume heteroscedastic measurement error.
#' @param bayes_factor logical; if TRUE, the Bayes factor is computed between the full model, and a model containing only the non-interaction surface.
#' @param robust logical: if TRUE, the model assumes a log-Pareto-tailed Normal (LPTN) distribution
#' @param rho numeric: the rho hyperparameter for the robust likelihood distribution
#' @param pcprior logical: if TRUE, the model uses a penalized complexity (PC) prior on the kernel hyperparameters of the Matérn kernel (type=3)
#' @param pcprior_hypers vector of size 4 giving hyperparameters for the PC prior on Matérn hyperparameters
#' @param lambda numeric; the parameter controls the residual noise observed in the heteroscedastic model when f = 0.
#' @param nu numeric; the nu parameter for the Matérn kernel. Must be one of (0.5, 1.5, 2.5)
#' @param method The method of estimation. Must be one of {`sampling`,`vb`} corresponding to full sampling, or variational Bayes.
#' @param control list; passed on to the stan sampler, e.g. for setting adapt_delta.
#' @param ... Arguments passed to \code{\link[rstan:sampling]{rstan::sampling}} or \code{\link[rstan:vb]{rstan::vb}} (e.g. iter, chains).
#'
#' @return An object of S3 class "\code{bayesynergy}", which is a list with the following entries
#' \tabular{ll}{
#' stanfit \tab An object of class \code{\link[rstan:stanmodel-class]{stanmodel}}, returned from the sampler. \cr
#' posterior_mean \tab A list containing the posterior means of model parameters. \cr
#' data \tab A list containing the original data used to fit the model. \cr
#' model \tab A list containing model specification for the model fit. \cr
#' returnCode \tab numeric; non-zero values indicate model was fit with errors or warnings \cr
#' LPML \tab numeric; The log pseudo-marginal likelihood of the fitted model. \cr
#' divergent \tab numeric; The number of divergent transitions from the Hamiltonian sampler. \cr
#' messages \tab character; any warnings are collected here for inspection. \cr
#' bayesfactor \tab numeric; The Bayes factor comparing the full model to a simpler model containing only the non-interaction assumption, p0. A large number indicating that an interaction effect is present. See the package vignette for more information.
#' }
#'
#'
#' @examples
#' \dontrun{
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
#' x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
#' fit <- bayesynergy(y_mat,x_mat)
#' }
#'
#' @importFrom utils modifyList
#' @importFrom bridgesampling bridge_sampler bf
#'
#' @export


bayesynergy <- function(y, x, type = 3, drug_names=NULL, experiment_ID = NULL, units = NULL,
                        lower_asymptotes = T, heteroscedastic = T, bayes_factor = F,
                        robust = F, rho = 0.90, pcprior = F, pcprior_hypers = c(1,0.1,1,0.2),
                         lambda = 0.005, nu = 1.5 , method = "sampling",
                        control = list(), ...){
  # Keep original data
  y.original = y
  x.original = x
  # Check that method is valid
  if (!(method %in% c("sampling","vb"))){
    stop("Method must be one of {'sampling','vb'}")
  }
  # Check that model specification is valid
  if (!(type%in%c(1,2,3,4))) {
    stop("Argument 'type' must be one of {1,2,3,4}!")
  }
  if (type == 2 & !(nu %in% c(0.5,1.5,2.5))){
    stop("Argument 'nu' must be one of {1/2,3/2,5/2}!")
  }
  if (type == 1 & bayes_factor){
    stop("Bayes factor calculation not available for type = 1")
  }
  if (type != 3 & pcprior){
    stop("PC prior only available for Matern covariance")
  }
  # Check that rho is correctly specified
  if (robust){
    if ((rho < (2*pnorm(1)-1)) | (rho > 1)){
      stop("rho must be between 0.6827 and 1")
    }
  }

  # Checking that data input is valid
  if (!is.numeric(x)){
    if (!(ncol(x)==2)){
      stop("Argument 'x' should be a two-column matrix of drug concentrations")
    }
  }
  if (nrow(as.matrix(y))!=nrow(x)){
    stop("Dimension mismatch! Arguments 'y' and 'x' should be equally sized.")
  }
  if (ncol(as.matrix(x))!=2){
    stop("Dimension mismatch! Argument 'x' should be a matrix with two columns containing drug concentrations.")
  }
  if (min(x[,1]) > 0 | min(x[,2]) > 0){
    stop("Missing monotherapy data! Make sure 'x' contains zero concentrations")
  }
  # Bayes factor requires sampling
  if (bayes_factor & method!="sampling"){
    stop("Computing the Bayes factor requires method = 'sampling')")
  }

  # GP models need adapt_delta > 0.9, so set that here
  control.default = list(adapt_delta = 0.8)
  if (type != 1){
    control.default = list(adapt_delta = 0.9)
  }
  control = modifyList(control.default,control)
  # Setting default names for drugs and experimentID
  if (is.null(drug_names)){
    drug_names = c("Drug1", "Drug2")
  }
  if (is.null(experiment_ID)){
    experiment_ID = "Experiment"
  }
  if (is.null(units)){
    units = c("conc.","conc.")
  }

  # Setting up data for the sampler
  unqX1 = log10(sort(unique(x[,1])))[-1] # Removing -Inf here
  unqX2 = log10(sort(unique(x[,2])))[-1] # Removing -Inf here

  # Need to find coordinates for the observed variables in this new coordinate system
  X1 = c(0,10^unqX1)
  X2 = c(0,10^unqX2)
  Xgrid = expand.grid(X1,X2)
  Xgrid = Xgrid[order(Xgrid[,"Var1"],Xgrid[,"Var2"]),]
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
    nmissing = (length(unqX1)+length(unqX2)+length(unqX1)*length(unqX2)+1)*nrep-length(ii_obs)
  }

  y = as.vector(y)
  y = y[which(!is.na(y))] # Removing missing

  # Create knots for the spline model
  dx1 = mean(diff(unqX1))
  dx2 = mean(diff(unqX2))
  n_knots1 = 5
  n_knots2 = 5
  degree = 3
  t1 = seq(min(unqX1)-0.5*dx1,max(unqX1)+0.5*dx1,length.out = n_knots1)
  t2 = seq(min(unqX2)-0.5*dx2,max(unqX2)+0.5*dx2,length.out = n_knots2)


  # Setting up data for Stan
  stan_data = list(n1=length(unqX1), n2=length(unqX2), x1=unqX1, x2=unqX2, nrep=nrep,
                   y=y, nmissing=nmissing, ii_obs = ii_obs, est_la = lower_asymptotes,
                   heteroscedastic = heteroscedastic, robust = robust, pcprior = pcprior,
                   lambda = lambda, pcprior_hypers = pcprior_hypers, rho = rho)
  stan_data_nointeraction = stan_data
  if (type == 1){ # Splines
    stan_data$n_knots1 = n_knots1
    stan_data$n_knots2 = n_knots2
    stan_data$t1 = t1
    stan_data$t2 = t2
    stan_data$degree = degree
  }
  else if (type == 2){ # GP w/  RBF kernel
    stan_data$kernel = 1
    stan_data$nu_matern = 0
    stan_data$est_alpha = 0
  }
  else if (type == 3){ # GP w/  Matérn kernel kernel
    stan_data$kernel = 2
    if (nu == 0.5){
      stan_data$nu_matern = 1
    } else if(nu == 1.5){
      stan_data$nu_matern = 2
    } else if (nu == 2.5){
      stan_data$nu_matern = 3
    }
    stan_data$est_alpha = 0
  }
  else if (type == 4){ # GP w/  RQ kernel
    stan_data$kernel = 3
    stan_data$nu_matern = 0
    stan_data$est_alpha = 1
  }

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
      if (type==1){
        fit = rstan::sampling(stanmodels$splines,stan_data, control = control, ...)
      }
      else {
        fit = rstan::sampling(stanmodels$gp_grid,stan_data, control = control, ...)
        if (bayes_factor){
          fitH0 = rstan::sampling(stanmodels$nointeraction,stan_data_nointeraction, control = control, ...)
        }
      }
    } else if (method == "vb"){
      if (type==1){
        fit = rstan::vb(stanmodels$splines,stan_data, ...)
      }
      else {
        fit = rstan::vb(stanmodels$gp_grid,stan_data, ...)
      }
    } else if (method == "opt"){
      if (type==1){
        fit = rstan::optimizing(stanmodels$splines,stan_data, as_vector = F, ...)
      }
      else {
        fit = rstan::optimizing(stanmodels$gp_grid,stan_data, as_vector = F, ...)
      }
    }
  )


  # Some calculations:
  posterior = rstan::extract(fit)
  n.save = length(posterior$lp__)
  LPML = -sum(log(apply(posterior$CPO,2,sum)/n.save))
  coef_names = names(posterior)
  # Remove those we don't care about
  coef_names = setdiff(coef_names,c("z","p0","p01","p02","Delta","wout","CPO","lp__"))

  # Surfaces
  # Mean response
  f = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  f[,1,1] = 1
  f[,2:(length(unqX2)+1),1] = posterior$p02
  f[,1,2:(length(unqX1)+1)] = posterior$p01
  f[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$p0+posterior$Delta
  f_mean = c()
  if (robust){
    f_mean = apply(f,c(2,3),median)
  } else{
    f_mean = apply(f,c(2,3),mean)
  }
  colnames(f_mean) = signif(c(0,10^unqX1),4)
  rownames(f_mean) = signif(c(0,10^unqX2),4)
  # Mean non-interaction
  p0 = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  p0[,1,1] = 1
  p0[,2:(length(unqX2)+1),1] = posterior$p02
  p0[,1,2:(length(unqX1)+1)] = posterior$p01
  p0[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$p0
  #p0_mean = apply(p0,c(2,3),mean)
  p0_mean = apply(p0,c(2,3),median)
  colnames(p0_mean) = signif(c(0,10^unqX1),4)
  rownames(p0_mean) = signif(c(0,10^unqX2),4)
  # Mean interaction
  Delta = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  Delta[,1,1] = 0
  Delta[,2:(length(unqX2)+1),1] = 0
  Delta[,1,2:(length(unqX1)+1)] = 0
  Delta[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$Delta
  #Delta_mean = apply(Delta,c(2,3),mean)
  Delta_mean = apply(Delta,c(2,3),median)
  colnames(Delta_mean) = signif(c(0,10^unqX1),4)
  rownames(Delta_mean) = signif(c(0,10^unqX2),4)

  # Pull out posterior mean and add some stuff
  posterior_mean = as.list(rstan::summary(fit,pars=coef_names)$summary[,'mean'])
  posterior_mean$f = f_mean
  posterior_mean$p0 = p0_mean
  posterior_mean$Delta = Delta_mean

  data = list(y = y.original, x = x.original, drug_names = drug_names, experiment_ID = experiment_ID, units = units, indices = ii_obs)
  model = list(type = type, lower_asymptotes = lower_asymptotes, method = method,
               bayes_factor = bayes_factor, heteroscedastic = heteroscedastic,
               robust = robust, pcprior = pcprior, lambda = lambda)

  # If matern kernel, add value of nu
  if (type==3){
    model$nu = nu
    model$pcprior_hypers = pcprior_hypers
  }

  # Inspect residuals
  residuals = y - as.vector(f_mean)[ii_obs]
  # Pull out pointwise observation noise
  point_sd = sqrt(posterior_mean$s^2*(f_mean[ii_obs]+lambda))
  # If any residuals are outside 3*sd, user is alerted
  if (sum(abs(residuals) > 3*point_sd) > 0){
    res_msg = paste("Largest residuals from posterior median is",signif(max(residuals),4),", which is more than three times the observation noise. This could indicate the presence of an outlier.")
    warning(res_msg,call. = F)
    returnCode = 1
    messages = c(messages,res_msg)
  }




  object = list(stanfit = fit, posterior_mean = posterior_mean,
                data = data, model = model, returnCode = returnCode,
                LPML = LPML)
  # Set some diagnostics
  object$divergent = NA
  if (returnCode){
    object$messages = messages
    if (method=="sampling"){
      object$divergent = sum(rstan::get_divergent_iterations(fit))
    }
  }

  # Calculate the Bayes Factor
  if (bayes_factor){
    message("Calculating the Bayes factor")
    full = bridge_sampler(fit,silent=T)
    no_interaction = bridge_sampler(fitH0,silent=T)
    bfactor = bf(full,no_interaction)
    object$bayesfactor = bfactor$bf
  }


  class(object) <- "bayesynergy"
  return(object)
}
