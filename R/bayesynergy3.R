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
    nmissing = (length(unqX1)+length(unqX2)+  length(unqX3) +
                  length(unqX1)*length(unqX2) + length(unqX1)*length(unqX3) +
                  length(unqX2)*length(unqX3) +
                  length(unqX1)*length(unqX2)*length(unqX3)+1)*nrep-length(ii_obs)
  }

  y = as.vector(y)
  y = y[which(!is.na(y))] # Removing missing


  # Setting up data for Stan
  stan_data = list(n1=length(unqX1), n2=length(unqX2), n3=length(unqX3),
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
  # posterior = rstan::extract(fit)
  # n.save = length(posterior$lp__)
  # LPML = -sum(log(apply(posterior$CPO,2,sum)/n.save))
  # coef_names = names(posterior)
  # # Remove those we don't care about
  # coef_names = setdiff(coef_names,c("z","p0","p01","p02","Delta","wout","CPO","lp__"))
  #
  # # Surfaces
  # # Mean response
  # f = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  # f[,1,1] = 1
  # f[,2:(length(unqX2)+1),1] = posterior$p02
  # f[,1,2:(length(unqX1)+1)] = posterior$p01
  # f[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$p0+posterior$Delta
  # f_mean = c()
  # if (robust){
  #   f_mean = apply(f,c(2,3),median)
  # } else{
  #   f_mean = apply(f,c(2,3),mean)
  # }
  # colnames(f_mean) = signif(c(0,10^unqX1),4)
  # rownames(f_mean) = signif(c(0,10^unqX2),4)
  # # Mean non-interaction
  # p0 = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  # p0[,1,1] = 1
  # p0[,2:(length(unqX2)+1),1] = posterior$p02
  # p0[,1,2:(length(unqX1)+1)] = posterior$p01
  # p0[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$p0
  # #p0_mean = apply(p0,c(2,3),mean)
  # p0_mean = apply(p0,c(2,3),median)
  # colnames(p0_mean) = signif(c(0,10^unqX1),4)
  # rownames(p0_mean) = signif(c(0,10^unqX2),4)
  # # Mean interaction
  # Delta = array(data=NA,c(n.save,length(unqX2)+1,length(unqX1)+1))
  # Delta[,1,1] = 0
  # Delta[,2:(length(unqX2)+1),1] = 0
  # Delta[,1,2:(length(unqX1)+1)] = 0
  # Delta[,2:(length(unqX2)+1),2:(length(unqX1)+1)] = posterior$Delta
  # #Delta_mean = apply(Delta,c(2,3),mean)
  # Delta_mean = apply(Delta,c(2,3),median)
  # colnames(Delta_mean) = signif(c(0,10^unqX1),4)
  # rownames(Delta_mean) = signif(c(0,10^unqX2),4)
  #
  # # Pull out posterior mean and add some stuff
  # posterior_mean = as.list(rstan::summary(fit,pars=coef_names)$summary[,'mean'])
  # posterior_mean$f = f_mean
  # posterior_mean$p0 = p0_mean
  # posterior_mean$Delta = Delta_mean
  #
  # data = list(y = y.original, x = x.original, drug_names = drug_names, experiment_ID = experiment_ID, units = units, indices = ii_obs)
  # model = list(type = type, lower_asymptotes = lower_asymptotes, method = method,
  #              bayes_factor = bayes_factor, heteroscedastic = heteroscedastic,
  #              robust = robust, pcprior = pcprior, lambda = lambda)
  #
  # # If matern kernel, add value of nu
  # if (type==3){
  #   model$nu = nu
  #   model$pcprior_hypers = pcprior_hypers
  # }
  #
  # # Inspect residuals
  # residuals = y - as.vector(f_mean)[ii_obs]
  # # Pull out pointwise observation noise
  # point_sd = sqrt(posterior_mean$s^2*(f_mean[ii_obs]+lambda))
  # # If any residuals are outside 3*sd, user is alerted
  # if (sum(abs(residuals) > 3*point_sd) > 0){
  #   res_msg = paste("Largest residuals from posterior median is",signif(max(residuals),4),", which is more than three times the observation noise. This could indicate the presence of an outlier.")
  #   warning(res_msg,call. = F)
  #   returnCode = 1
  #   messages = c(messages,res_msg)
  # }
  #
  #
  #
  #
  # object = list(stanfit = fit, posterior_mean = posterior_mean,
  #               data = data, model = model, returnCode = returnCode,
  #               LPML = LPML)
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
  # class(object) <- "bayesynergy"
  return(fit)
}
