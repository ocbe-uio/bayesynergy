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
#' @param log10_conc logical; if TRUE concentrations are assumed given on the log10 scale.
#' @param lower_asymptotes logical; if TRUE the model will estimate the lower asymptotes of monotherapy curves.
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
#' }
#' 
#' 
#' @examples 
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
#' x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
#' fit <- bayesynergy(y_mat,x_mat)
#' 
#' @export


bayesynergy <- function(y, x, type = 3, drug_names=NULL, experiment_ID = NULL, log10_conc = FALSE, 
                        lower_asymptotes = T, nu = 1.5 , method = "sampling",
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
  if (log10_conc){
    x = 10^x
  }
  if (min(x[,1]) > 0 | min(x[,2]) > 0){
    stop("Missing monotherapy data! Make sure 'x' contains zero concentrations (-Inf if given on log10 scale).")
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
  
  # Setting up data for the sampler
  unqX1 = log10(sort(unique(x[,1])))[-1] # Removing -Inf here
  unqX2 = log10(sort(unique(x[,2])))[-1] # Removing -Inf here
  
  # Need to find coordinates for the observed variables in this new coordinate system
  X1 = c(0,10^unqX1)
  X2 = c(0,10^unqX2)
  Xgrid = expand.grid(X1,X2)
  Xgrid = Xgrid[order(Xgrid["Var1"],Xgrid["Var2"]),]
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
       pij=y, nmissing=nmissing, ii_obs = ii_obs, est_la = lower_asymptotes)
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
        fit = rstan::sampling(stanmodels$splines,stan_data, pars=c("z"),include=F, control = control, ...)
      }
      else {
        fit = rstan::sampling(stanmodels$gp_grid,stan_data, pars=c("z"),include=F, control = control, ...)
      }
    } else if (method == "vb"){
      if (type==1){
        fit = rstan::vb(stanmodels$splines,stan_data, pars=c("z"),include=F, ...)
      }
      else {
        fit = rstan::vb(stanmodels$gp_grid,stan_data, pars=c("z"),include=F, ...)
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
  coef_names = setdiff(coef_names,c("pij_0","pij_01","pij_02","Delta_ij","CPO","lp__"))

  # Surfaces
  # Mean response
  p_ij = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  p_ij[,1,1] = 1
  p_ij[,2:(length(unqX2)+1),1] = posterior$pij_02
  p_ij[,1,2:(length(unqX1)+1)] = posterior$pij_01
  p_ij[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$pij_0+posterior$Delta_ij
  p_ij_mean = apply(p_ij,c(2,3),mean)
  colnames(p_ij_mean) = round(c(0,10^unqX1),4)
  rownames(p_ij_mean) = round(c(0,10^unqX2),4)
  # Mean non-interaction
  p_0 = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  p_0[,1,1] = 1
  p_0[,2:(length(unqX2)+1),1] = posterior$pij_02
  p_0[,1,2:(length(unqX1)+1)] = posterior$pij_01
  p_0[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$pij_0
  p_0_mean = apply(p_0,c(2,3),mean)
  colnames(p_0_mean) = round(c(0,10^unqX1),4)
  rownames(p_0_mean) = round(c(0,10^unqX2),4)
  # Mean interaction
  Delta = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  Delta[,1,1] = 0
  Delta[,2:(length(unqX2)+1),1] = 0
  Delta[,1,2:(length(unqX1)+1)] = 0
  Delta[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$Delta_ij
  Delta_mean = apply(Delta,c(2,3),mean)
  colnames(Delta_mean) = round(c(0,10^unqX1),4)
  rownames(Delta_mean) = round(c(0,10^unqX2),4)
  
  
  posterior_mean = as.list(rstan::summary(fit,pars=coef_names)$summary[,'mean'])
  posterior_mean$p_ij = p_ij_mean
  posterior_mean$p_0 = p_0_mean
  posterior_mean$Delta = Delta_mean
  
  data = list(y = y.original, x = x.original, drug_names = drug_names, experiment_ID = experiment_ID)
  model = list(type = type, lower_asymptotes = lower_asymptotes, method = method)
  
  # If matern kernel, add value of nu
  if (type==3){
    model$nu = nu
  }
  
  object = list(stanfit = fit, posterior_mean = posterior_mean, 
                data = data, model = model, returnCode = returnCode,
                LPML = LPML)
  if (returnCode){
    object$messages = messages
  }
  
  class(object) <- "bayesynergy"
  return(object)
}
  