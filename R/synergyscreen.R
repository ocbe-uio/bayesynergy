#' Function for fitting high-througput drug combination screens with parallel processing
#' 
#' @description 
#' The function \code{synergyscreen} allows the fitting of high-throughput drug combination screens through parallel processing
#' 
#'    
#' @param experiments A list of experiments obtained from a high-throughput screen. See *Details* for more information on the structure of each element.
#' @param metric A vector of metrics of interest returned to the user. Must be a vector of named variables from the \code{posterior_mean} list of the S3 \code{bayesynergy} object.
#' @param return_samples logical; if TRUE, the function returns posterior samples of each metric.
#' @param save_raw logical; if TRUE, the raw bayesynergy object is saved for each individual experiment.
#' @param save_plots  logical; if TRUE, plots for each individual experiment is saved.
#' @param path string; path for saving output and plot for each individual experiment.
#' @param parallel logical; if TRUE, parallel processing is utilized to run the screen.
#' @param max_cores integer; the maximum number of cores to utilize for the parallel processing.
#' @param plot_params list; parameters to be passed to the plotting function. See \link{plot.bayesynergy} for details.
#' @param bayesynergy_params list; parameters to be passed to the bayesynergy function. See \link{bayesynergy} for details.
#'  
#' 
#' @details
#' The elements of \strong{experiments} must themselves be lists with the following elements
#' \tabular{ll}{
#' y: \tab matrix of viability scores \cr
#' x: \tab matrix of concentrations \cr
#' drug_names: \tab vector of drug names in the experiment \cr
#' experiment_ID: \tab string denoting the unique experiment ID, e.g. cell line name.
#' }
#' 
#' @return A list of equal length as \code{experiments}, each element containing
#' \tabular{ll}{
#' summaries \tab summary statistics for variables defined in \code{metric}. \cr
#' drug_names \tab names of the drugs utilized for the experiment. \cr
#' experiment_ID \tab identifier of experiment, typically name of cell Line. \cr
#' samples \tab if requested, posterior samples of variables defined in \code{metric}.
#' }
#' 
#' 
#'
#' @examples
#' \dontrun{
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' experiment1 = list(y = mathews_DLBCL[[1]][[1]], 
#' x = mathews_DLBCL[[1]][[2]], 
#' drug_names = c("ispinesib","ibrutinib"))
#' experiment2 = list(y = mathews_DLBCL[[2]][[1]], 
#' x = mathews_DLBCL[[2]][[2]], 
#' drug_names = c("canertinib","ibrutinib"))
#' experiments = list(experiment1,experiment2)
#' fit <- synergyscreen(experiments)
#' }
#'
#' @export
#'
#' @import foreach doParallel parallel doSNOW
#' @importFrom utils setTxtProgressBar txtProgressBar


synergyscreen = function(experiments, metric = c("rVUS_syn","rVUS_ant"), return_samples = F,
                         save_raw = T, save_plots = T, path = NULL, parallel=T, max_cores=NULL,
                         plot_params = list(), bayesynergy_params = list()){
  # Check that path is not null, and if so, set it to work directory
  if (is.null(path)){
    path <- getwd()
  }
  # Check that path exists, if not create it
  if (!dir.exists(path)){
    dir.create(path)
  }
  # Put dirmark on file path
  path = Sys.glob(path,dirmark = T)
  # Create container for results
  results = c()
  if (parallel){
    # Setup cluster
    # Some diversion here for CRAN checks
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      max_cores <- 2L
    } else {
      # use all cores in devtools::test()
      max_cores <- min(length(experiments),min(max_cores, parallel::detectCores() - 1))
    }
    print(paste("Using",max_cores,"cores for parallel run on list of size",length(experiments)))
    print(paste("Saving output at:",path))
    cl = parallel::makeCluster(max_cores,setup_strategy = "sequential")
    registerDoSNOW(cl)
    # Setting up progress bar
    pb <- txtProgressBar(max = length(experiments), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    # Running the experiments, in parallel
    results = foreach (ee=experiments, .combine=list,
                       .maxcombine = length(experiments),
                       .multicombine = T,
                       .options.snow = opts,
                       .packages = "bayesynergy") %dopar% {
                         # Do the fitting here
                         data <- ee
                         # If drug_names or experiment_ID not given
                         if (!("drug_names" %in% names(data))){
                           data$drug_names = c("DrugA","DrugB")
                         }
                         if (!("experiment_ID" %in% names(data))){
                           data$experiment_ID = "Experiment"
                         }
                         fit <- do.call(bayesynergy,c(data,bayesynergy_params))
                         # Saving raw
                         if (save_raw){
                           save(fit,file=paste0(path,paste(data$experiment_ID,data$drug_names[1],data$drug_names[2],"raw",sep="_")))
                         }
                         # Saving plots
                         if (save_plots){
                           suppressMessages(do.call(plot,c(list(x=fit,save_plot = T, path=path), plot_params)))
                         }
                         
                         # Posterior
                         toReturn = rstan::summary(fit$stanfit)$summary[metric,]
                         samples = rstan::extract(fit$stanfit,pars=metric)
                         if (return_samples){
                           list(summaries = toReturn, drug_names = data$drug_names, experiment_ID=data$experiment_ID,samples=samples)
                         } else {
                           list(summaries = toReturn, drug_names = data$drug_names, experiment_ID=data$experiment_ID)
                         }
                       }
    # Tidying up the progress bar
    close(pb)
    # Stopping cluster
    stopCluster(cl)
  } else { # Simple for loop
    for (i in 1:length(experiments)){
      data <- experiments[[i]]
      # If drug_names or experiment_ID not given
      if (!("drug_names" %in% names(data))){
        data$drug_names = c("DrugA","DrugB")
      }
      if (!("experiment_ID" %in% names(data))){
        data$experiment_ID = "Experiment"
      }
      fit <- do.call(bayesynergy,c(data,bayesynergy_params))
      # Saving raw
      if (save_raw){
        save(fit,file=paste0(path,paste(data$experiment_ID,data$drug_names[1],data$drug_names[2],"raw",sep="_")))
      }
      # Saving plots
      if (save_plots){
        suppressMessages(do.call(plot,c(list(x=fit,save_plot = T, path=path), plot_params)))
      }
      
      # Posterior
      toReturn = rstan::summary(fit$stanfit)$summary[metric,]
      samples = rstan::extract(fit$stanfit,pars=metric)
      if (return_samples){
        results[[i]] = list(summaries = toReturn, drug_names = data$drug_names, experiment_ID=data$experiment_ID,samples=stats)
      } else {results[[i]] = list(summaries = toReturn, drug_names = data$drug_names, experiment_ID=data$experiment_ID)}
    }
  }
  
  # Returning
  return(results)
  
}

