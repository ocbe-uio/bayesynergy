#' Function for fitting high-througput drug combination screens with parallel processing
#' 
#' @description 
#' The function \code{synergyscreen} allows the fitting of high-throughput drug combination screens through parallel processing
#' 
#'  
#'    
#' @param experiments A list of experiments obtained from a high-throughput screen. See *Details* for more information on the structure of each element
#' @param metric The metric returned to the user. Must be a named variable from the \code{Summary_Ouput} of the S3 BayeSyneRgy object.
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
#' y: \tab matrix of concentrations \cr
#' drug_names: \tab vector of drug names in the experiment \cr
#' experiment_ID: \tab string denoting the unique experiment ID, e.g. cell line name.
#' }
#' 
#' 
#' 
#'
#' @examples
#' \dontrun{
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' experiment1 = list(y = mathews_DLBCL[[1]][[1]], x = mathews_DLBCL[[1]][[2]], drug_names = c("ispinesib","ibrutinib"),experiment_ID = "experiment1")
#' experiment2 = list(y = mathews_DLBCL[[2]][[1]], x = mathews_DLBCL[[2]][[2]], drug_names = c("canertinib","ibrutinib"),experiment_ID = "experiment2")
#' experiments = list(experiment1,experiment2)
#' fit <- synergyscreen(experiments)
#' }
#'
#' @export
#'
#' @import foreach doParallel parallel doSNOW

synergyscreen = function(experiments, metric = c("rVUS_syn","rVUS_ant"), save_raw = T, save_plots = T, path = NULL, parallel=T, max_cores=NULL,
                         plot_params = list(), bayesynergy_params = list()){
  # Check that path is not null, and if so, set it to work directory
  if (is.null(path)){
    path <- getwd()
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
                         posterior = rstan::extract(fit$stanfit)
                         toReturn = matrix(NA,ncol=5,nrow=length(metric))
                         rownames(toReturn) = metric
                         colnames(toReturn) = c("mean","sd","2.5%","50%","97.5%")
                         stats = list()
                         for (m in  1:length(metric)){
                           statistic = posterior[[metric[m]]]
                           summaries = c(mean(statistic),sd(statistic),quantile(statistic,probs=c(0.025,0.5,0.975)))
                           stats[[m]] = list(name=metric[m],samples = statistic)
                           toReturn[m,] = summaries
                         }
                         list(summaries = toReturn, drug_names = data$drug_names, experiment_ID=data$experiment_ID,samples=stats)
                       }
    # Tidying up the progress bar
    close(pb)
    # Stopping cluster
    stopCluster(cl)
  } else { # Simple for loop
    for (i in 1:length(experiments)){
      data <- experiments[[i]]
      fit <- do.call(bayesynergy,c(data,bayesynergy_params))
      # Saving raw
      if (save_raw){
        save(fit,file=paste0(path,paste(data$experiment_ID,data$drug_names[1],data$drug_names[2],"raw",sep="_")))
      }
      # Saving plots
      if (save_plots){
        suppressMessages(do.call(plot,c(fit,save_plot = T, path=path, plot_params)))
      }
      
      # Posterior 
      posterior = rstan::extract(fit$stanfit)
      toReturn = matrix(NA,ncol=5,nrow=length(metric))
      rownames(toReturn) = metric
      colnames(toReturn) = c("mean","sd","2.5%","50%","97.5%")
      stats = list()
      for (m in  1:length(metric)){
        statistic = posterior[[metric[m]]]
        summaries = c(mean(statistic),sd(statistic),quantile(statistic,probs=c(0.025,0.5,0.975)))
        stats[[m]] = list(name=metric[m],samples = statistic)
        toReturn[m,] = summaries
      }
      results[[i]] = list(summaries = toReturn, drug_names = data$drug_names, experiment_ID=data$experiment_ID,samples=stats)
    }
  }
  
  # Returning
  return(results)
  
}

