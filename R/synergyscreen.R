#' Function for fitting high-througput drug combination screens with parallel processing
#' 
#' @description 
#' The function \code{synergyscreen} allows the fitting of high-throughput drug combination screens through parallel processing
#' 
#'    
#' @param experiments A list of experiments obtained from a high-throughput screen. See *Details* for more information on the structure of each element.
#' @param return_samples logical; if TRUE, the function returns the full fitted \code{\link[bayesynergy]{bayesynergy}} object.
#' @param save_raw logical; if TRUE, the raw bayesynergy object is saved for each individual experiment.
#' @param save_plots  logical; if TRUE, plots for each individual experiment is saved.
#' @param path string; path for saving output and plot for each individual experiment.
#' @param parallel logical; if TRUE, parallel processing is utilized to run the screen.
#' @param max_cores integer; the maximum number of cores to utilize for the parallel processing.
#' @param max_retries intereger; the maximum number of retries utilized in model fit.
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
#' summaries \tab posterior summary statistics for variables of the model. \cr
#' drug_names \tab names of the drugs utilized for the experiment. \cr
#' experiment_ID \tab identifier of experiment, typically name of cell Line. \cr
#' fit \tab if requested, the fitted \code{\link[bayesynergy]{bayesynergy}} object.
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


synergyscreen = function(experiments, return_samples = F,
                         save_raw = T, save_plots = T, path = NULL, parallel=T, max_cores=NULL,
                         max_retries = 3,
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
  # Telling the user where things are being saved
  if (save_raw | save_plots){
    print(paste("Saving output at:",path))
  }
  
  # We need to check that bayesyn_params contain some parameters
  if (!("control" %in% names(bayesynergy_params))){
    bayesynergy_params$control = list()
  } 
  if (!("method" %in% names(bayesynergy_params))){
    bayesynergy_params$method = "sampling"
  }
  
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
                       .packages = "bayesynergy",
                       .verbose=F,
                       .errorhandling = "pass") %dopar% {
                         # Do the fitting here
                         data <- ee
                         # If drug_names or experiment_ID not given
                         if (!("drug_names" %in% names(data))){
                           data$drug_names = c("DrugA","DrugB")
                         }
                         if (!("experiment_ID" %in% names(data))){
                           data$experiment_ID = "Experiment"
                         }
                         # Wrap this in a try-catch
                         fit = c()
                         suppressWarnings({
                           tryCatch(
                             {fit = do.call(bayesynergy,c(data,bayesynergy_params))},
                             error = function(e){fit$returnCode <<- 2}
                           )
                           retry = T
                           retries = 0
                           while (retry){
                             # Retry again with stricter settings
                             if (bayesynergy_params$method == "sampling"){ # For sampler, we set adapt delta to high value and hope for the best
                               if (fit$returnCode != 0){
                                 bayesynergy_params$control <- modifyList(bayesynergy_params$control,list(adapt_delta=0.99))
                                 tryCatch(
                                   {fit = do.call(bayesynergy,c(data,bayesynergy_params))},
                                   error = function(e){fit$returnCode <<- 2}
                                 )
                               }
                             } else if (bayesynergy_params$method == "vb"){ 
                               # For VB, we simply need to try again if the previous one gave an error
                               # or clearly has a terrible fit
                               if ((fit$returnCode == 2) || (fit$posterior_mean$s > 1)){
                                 tryCatch(
                                   {fit = do.call(bayesynergy,c(data,bayesynergy_params))},
                                   error = function(e){fit$returnCode <<- 2}
                                 )
                               }
                             }
                             retries = retries + 1
                             if (retries > max_retries){retry = F}
                           }
                           
                           
                           if (fit$returnCode != 2){
                             # Saving raw
                             if (save_raw){
                               save(fit,file=paste0(path,paste(data$experiment_ID,data$drug_names[1],data$drug_names[2],"raw",sep="_")))
                             }
                             # Saving plots
                             if (save_plots){
                               suppressMessages(do.call(plot,c(list(x=fit,save_plot = T, path=path), plot_params)))
                             }
                             
                             # Create some summaries
                             summaryStats = rstan::summary(fit$stanfit)$summary
                             synMetrics <- data.frame(
                               `Experiment ID` = data$experiment_ID,
                               `Drug A` = data$drug_names[1],
                               `Drug B` = data$drug_names[2],
                               # Extract the EC50s for each drug
                               `EC50 (Drug A)` = summaryStats["ec50_1","mean"],
                               `EC50 (Drug B)` = summaryStats["ec50_2","mean"],
                               
                               # Calculate a standardized synergy score from the mean of the rVUS
                               `Synergy Score` = (summaryStats["log_rVUS_syn","mean"]/summaryStats["log_rVUS_syn","sd"]),
                               
                               # Calculating additional statistics and parameters
                               `Mean (syn)` = summaryStats["rVUS_syn", "mean"],
                               `SEM (syn)` = summaryStats["rVUS_syn", "se_mean"],
                               `SD (syn)` = summaryStats["rVUS_syn", "sd"],
                               `Mean/SD (syn)` = summaryStats["rVUS_syn", "mean"] / summaryStats["rVUS_syn", "sd"],
                               `97.5% (syn)` = summaryStats["rVUS_syn", "97.5%"],
                               
                               `Antagonism Score` = (summaryStats["log_rVUS_ant","mean"]/summaryStats["log_rVUS_ant","sd"]),
                               `Mean (ant)` = summaryStats["rVUS_ant", "mean"],
                               `SEM (ant)` = summaryStats["rVUS_ant", "se_mean"],
                               `SD (ant)` = summaryStats["rVUS_ant", "sd"],
                               `Mean/SD (ant)` = summaryStats["rVUS_ant", "mean"] / summaryStats["rVUS_ant", "sd"],
                               `97.5% (ant)` = summaryStats["rVUS_ant", "97.5%"],
                               `s` = summaryStats["s", "mean"],
                               `returnCode` = fit$returnCode,
                               `divergentTransitions` = sum(fit$divergent),
                               
                               check.names = FALSE, stringsAsFactors = FALSE)
                           } else {
                             # Create some summaries
                             summaryStats = rstan::summary(fit$stanfit)$summary
                             synMetrics <- data.frame(
                               `Experiment ID` = data$experiment_ID,
                               `Drug A` = data$drug_names[1],
                               `Drug B` = data$drug_names[2],
                               # Extract the EC50s for each drug
                               `EC50 (Drug A)` = NA,
                               `EC50 (Drug B)` = NA,
                               
                               # Calculate a standardized synergy score from the mean of the rVUS
                               `Synergy Score` = NA,
                               
                               # Calculating additional statistics and parameters
                               `Mean (syn)` = NA,
                               `SEM (syn)` = NA,
                               `SD (syn)` = NA,
                               `Mean/SD (syn)` = NA,
                               `97.5% (syn)` = NA,
                               
                               `Antagonism Score` = NA,
                               `Mean (ant)` = NA,
                               `SEM (ant)` = NA,
                               `SD (ant)` = NA,
                               `Mean/SD (ant)` = NA,
                               `97.5% (ant)` = NA,
                               `s` = NA,
                               `returnCode` = fit$returnCode,
                               `divergentTransitions` = NA,
                               
                               check.names = FALSE, stringsAsFactors = FALSE)
                             
                           }
                         })
                         if (return_samples & synMetrics[,"returnCode"] != 2){
                           list(summaries = synMetrics, drug_names = data$drug_names, experiment_ID=data$experiment_ID,fit=fit)
                         } else {
                           list(summaries = synMetrics, drug_names = data$drug_names, experiment_ID=data$experiment_ID)
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
      # Wrap this in a try-catch
      fit = c()
      suppressWarnings({
        tryCatch(
          {fit = do.call(bayesynergy,c(data,bayesynergy_params))},
          error = function(e){fit$returnCode <<- 2}
        )
        retry = T
        retries = 0
        while (retry){
          # Retry again with stricter settings
          if (bayesynergy_params$method == "sampling"){ # For sampler, we set adapt delta to high value and hope for the best
            if (fit$returnCode != 0){
              bayesynergy_params$control <- modifyList(bayesynergy_params$control,list(adapt_delta=0.99))
              tryCatch(
                {fit = do.call(bayesynergy,c(data,bayesynergy_params))},
                error = function(e){fit$returnCode <<- 2}
              )
            }
          } else if (bayesynergy_params$method == "vb"){ 
            # For VB, we simply need to try again if the previous one gave an error
            # or clearly has a terrible fit
            if ((fit$returnCode == 2) || (fit$posterior_mean$s > 1)){
              tryCatch(
                {fit = do.call(bayesynergy,c(data,bayesynergy_params))},
                error = function(e){fit$returnCode <<- 2}
              )
            }
          }
          retries = retries + 1
          if (retries > max_retries){retry = F}
        }
        
        
        if (fit$returnCode != 2){
          # Saving raw
          if (save_raw){
            save(fit,file=paste0(path,paste(data$experiment_ID,data$drug_names[1],data$drug_names[2],"raw",sep="_")))
          }
          # Saving plots
          if (save_plots){
            suppressMessages(do.call(plot,c(list(x=fit,save_plot = T, path=path), plot_params)))
          }
          
          # Create some summaries
          summaryStats = rstan::summary(fit$stanfit)$summary
          synMetrics <- data.frame(
            `Experiment ID` = data$experiment_ID,
            `Drug A` = data$drug_names[1],
            `Drug B` = data$drug_names[2],
            # Extract the EC50s for each drug
            `EC50 (Drug A)` = summaryStats["ec50_1","mean"],
            `EC50 (Drug B)` = summaryStats["ec50_2","mean"],
            
            # Calculate a standardized synergy score from the mean of the rVUS
            `Synergy Score` = (summaryStats["log_rVUS_syn","mean"]/summaryStats["log_rVUS_syn","sd"]),
            
            # Calculating additional statistics and parameters
            `Mean (syn)` = summaryStats["rVUS_syn", "mean"],
            `SEM (syn)` = summaryStats["rVUS_syn", "se_mean"],
            `SD (syn)` = summaryStats["rVUS_syn", "sd"],
            `Mean/SD (syn)` = summaryStats["rVUS_syn", "mean"] / summaryStats["rVUS_syn", "sd"],
            `97.5% (syn)` = summaryStats["rVUS_syn", "97.5%"],
            
            `Antagonism Score` = (summaryStats["log_rVUS_ant","mean"]/summaryStats["log_rVUS_ant","sd"]),
            `Mean (ant)` = summaryStats["rVUS_ant", "mean"],
            `SEM (ant)` = summaryStats["rVUS_ant", "se_mean"],
            `SD (ant)` = summaryStats["rVUS_ant", "sd"],
            `Mean/SD (ant)` = summaryStats["rVUS_ant", "mean"] / summaryStats["rVUS_ant", "sd"],
            `97.5% (ant)` = summaryStats["rVUS_ant", "97.5%"],
            `s` = summaryStats["s", "mean"],
            `returnCode` = fit$returnCode,
            `divergentTransitions` = sum(fit$divergent),
            
            check.names = FALSE, stringsAsFactors = FALSE)
        } else {
          # Create some summaries
          summaryStats = rstan::summary(fit$stanfit)$summary
          synMetrics <- data.frame(
            `Experiment ID` = data$experiment_ID,
            `Drug A` = data$drug_names[1],
            `Drug B` = data$drug_names[2],
            # Extract the EC50s for each drug
            `EC50 (Drug A)` = NA,
            `EC50 (Drug B)` = NA,
            
            # Calculate a standardized synergy score from the mean of the rVUS
            `Synergy Score` = NA,
            
            # Calculating additional statistics and parameters
            `Mean (syn)` = NA,
            `SEM (syn)` = NA,
            `SD (syn)` = NA,
            `Mean/SD (syn)` = NA,
            `97.5% (syn)` = NA,
            
            `Antagonism Score` = NA,
            `Mean (ant)` = NA,
            `SEM (ant)` = NA,
            `SD (ant)` = NA,
            `Mean/SD (ant)` = NA,
            `97.5% (ant)` = NA,
            `s` = NA,
            `returnCode` = fit$returnCode,
            `divergentTransitions` = NA,
            
            check.names = FALSE, stringsAsFactors = FALSE)
          
        }
      })
    
      if (return_samples){
        results[[i]] = list(summaries = synMetrics, drug_names = data$drug_names, experiment_ID=data$experiment_ID,fit=fit)
      } else {results[[i]] = list(summaries = synMetrics, drug_names = data$drug_names, experiment_ID=data$experiment_ID)}
    }
  }
  
  # Combine and reorder some stuff
  toReturn = list()
  screenSummary = do.call(rbind,lapply(results,function(x) x$summaries))
  toReturn$screenSummary = screenSummary
  if (return_samples){
    screenSamples = lapply(results, function(x) x$fit)
    toReturn$screenSamples = screenSamples
  }
  
  # Do some checks to see if things converges okay
  propFailed = mean((toReturn$screenSummary[,"returnCode"]>0))
  if (propFailed > 0){
    if (bayesynergy_params$method == "sampling"){
      warning(paste0(round(100*propFailed,digits=2),"% of experiments returned a warning or error -- check these."))
    }
    if (bayesynergy_params$method == "vb"){
      warning(paste0(round(100*propFailed,digits=2),"% of experiments returned a warning or error. NOTE: For method = 'vb' this is to be expected"))
    }
  }
  if (nrow(toReturn$screenSummary) != length(experiments)){
    diff = length(experiments) - nrow(toReturn$screenSummary)
    warning(paste(diff," experiments failed to process"))
  }

  # Returning
  class(toReturn) <- "synergyscreen"
  return(toReturn)
  
}

