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
#' @param max_retries integer; the maximum number of retries utilized in model fit.
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
#' units vector of size 2; concentration units for the drugs, e.g. c("\eqn{\mu}M","\eqn{\mu}M")
#' }
#'
#' @return A list containing two elements
#' \tabular{ll}{
#' screenSummary \tab data frame, posterior summary statistics for each experiment. \cr
#' failed \tab A list containing elements experiments that failed to process. \cr
#' screenSamples \tab if requested, a list containing the fitted \code{\link[bayesynergy]{bayesynergy}} object for each experiment.
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

  # Default max cores
  if (is.null(max_cores)){
    max_cores <- min(length(experiments),min(max_cores, parallel::detectCores() - 1))
  } else { # Overwrite on user input if allowed
    max_cores = min(max_cores, parallel::detectCores())
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
  # Set default method to vb for computational speed
  if (!("method" %in% names(bayesynergy_params))){
    bayesynergy_params$method = "vb"
  }
  # Finally, if we are re-running failed experiments
  rerun = F
  oldrun = c()
  if (class(experiments)=="synergyscreen"){
    if ("failed" %in% names(experiments)){
      oldrun = experiments$screenSummary
      experiments = experiments$failed
      rerun = T
    } else (stop("No failed experiments found"))
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
    }

    print(paste("Using",max_cores,"cores for parallel run on list of size",length(experiments)))
    cl = parallel::makeCluster(max_cores,setup_strategy = "sequential")
    registerDoSNOW(cl)
    # Setting up progress bar
    pb <- txtProgressBar(max = length(experiments), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    # Running the experiments, in parallel
    results = foreach (kk = 1:length(experiments), .combine=list,
                       .maxcombine = length(experiments),
                       .multicombine = T,
                       .options.snow = opts,
                       .packages = "bayesynergy",
                       .verbose=F,
                       .errorhandling = "pass") %dopar% {
                         # Do the fitting here
                         data <- experiments[[kk]]
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
                               suppressMessages(do.call(plot,c(list(x=fit,save_plots = T, path=path), plot_params)))
                             }

                             # Create some summaries
                             summaryStats = rstan::summary(fit$stanfit)$summary
                             synMetrics <- data.frame(
                               `listID` = kk,
                               `Experiment ID` = data$experiment_ID,
                               `Drug A` = data$drug_names[1],
                               `Drug B` = data$drug_names[2],
                               # Include some summary stats for the experiment
                               `EC50 (Drug A)` = summaryStats["ec50_1","mean"],
                               `EC50 (Drug B)` = summaryStats["ec50_2","mean"],
                               `DSS (Drug A)` = summaryStats["dss_1","mean"],
                               `DSS (Drug B)` = summaryStats["dss_2","mean"],
                               `Synergy (mean)` = summaryStats["VUS_syn","mean"],
                               `Synergy (sd)` = summaryStats["VUS_syn","sd"],
                               `Antagonism (mean)` = summaryStats["VUS_ant","mean"],
                               `Antagonism (sd)` = summaryStats["VUS_ant","sd"],
                               `Efficacy (mean)` = summaryStats["rVUS_f","mean"],
                               `Efficacy (sd)` = summaryStats["rVUS_f","sd"],
                               `Non-interaction Efficacy (mean)` = summaryStats["rVUS_p0","mean"],
                               `Non-interaction Efficacy (sd)` = summaryStats["rVUS_p0","sd"],
                               `Interaction (mean)` = summaryStats["VUS_Delta","mean"],
                               `Interaction (sd)` = summaryStats["VUS_Delta","sd"],

                               # Calculate standardized scores for comparisons
                               `Synergy Score` = summaryStats["VUS_syn","mean"]/summaryStats["VUS_syn","sd"],
                               `Antagonism Score` = summaryStats["VUS_ant","mean"]/summaryStats["VUS_ant","sd"],
                               # Finally some information about model fit
                               `s` = summaryStats["s", "mean"],
                               `returnCode` = fit$returnCode,
                               `divergentTransitions` = sum(fit$divergent),

                               check.names = FALSE, stringsAsFactors = FALSE)
                             if (fit$model$bayes_factor){
                               synMetrics$`Bayes Factor` = fit$bayesfactor
                             }
                           } else {
                             # Create some summaries
                             summaryStats = rstan::summary(fit$stanfit)$summary
                             synMetrics <- data.frame(
                               `listID` = kk,
                               `Experiment ID` = data$experiment_ID,
                               `Drug A` = data$drug_names[1],
                               `Drug B` = data$drug_names[2],
                               # Include some summary stats for the experiment
                               `EC50 (Drug A)` = NA,
                               `EC50 (Drug B)` = NA,
                               `DSS (Drug A)` = NA,
                               `DSS (Drug B)` = NA,
                               `Synergy (mean)` = NA,
                               `Synergy (sd)` = NA,
                               `Antagonism (mean)` = NA,
                               `Antagonism (sd)` = NA,
                               `Efficacy (mean)` = NA,
                               `Efficacy (sd)` = NA,
                               `Non-interaction Efficacy (mean)` = NA,
                               `Non-interaction Efficacy (sd)` = NA,
                               `Interaction (mean)` = NA,
                               `Interaction (sd)` = NA,

                               # Calculate standardized scores for comparisons
                               `Synergy Score` = NA,
                               `Antagonism Score` = NA,
                               # Finally some information about model fit
                               `s` = NA,
                               `returnCode` = 2,
                               `divergentTransitions` = NA,

                               check.names = FALSE, stringsAsFactors = FALSE)
                             if (fit$model$bayes_factor){
                               synMetrics$`Bayes Factor` = NA
                             }
                           }
                         })

                         # Garbage collection
                         dir <- tempdir()
                         file.names <- list.files(path = dir, pattern="*.csv")
                         unlink(paste0(dir,"/",file.names))

                         # Return
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
            suppressMessages(do.call(plot,c(list(x=fit,save_plots = T, path=path), plot_params)))
          }

          # Create some summaries
          summaryStats = rstan::summary(fit$stanfit)$summary
          synMetrics <- data.frame(
            `listID` = i,
            `Experiment ID` = data$experiment_ID,
            `Drug A` = data$drug_names[1],
            `Drug B` = data$drug_names[2],
            # Include some summary stats for the experiment
            `EC50 (Drug A)` = summaryStats["ec50_1","mean"],
            `EC50 (Drug B)` = summaryStats["ec50_2","mean"],
            `DSS (Drug A)` = summaryStats["dss_1","mean"],
            `DSS (Drug B)` = summaryStats["dss_2","mean"],
            `Synergy (mean)` = summaryStats["VUS_syn","mean"],
            `Synergy (sd)` = summaryStats["VUS_syn","sd"],
            `Antagonism (mean)` = summaryStats["VUS_ant","mean"],
            `Antagonism (sd)` = summaryStats["VUS_ant","sd"],
            `Efficacy (mean)` = summaryStats["rVUS_f","mean"],
            `Efficacy (sd)` = summaryStats["rVUS_f","sd"],
            `Non-interaction Efficacy (mean)` = summaryStats["rVUS_p0","mean"],
            `Non-interaction Efficacy (sd)` = summaryStats["rVUS_p0","sd"],
            `Interaction (mean)` = summaryStats["VUS_Delta","mean"],
            `Interaction (sd)` = summaryStats["VUS_Delta","sd"],

            # Calculate standardized scores for comparisons
            `Synergy Score` = summaryStats["VUS_syn","mean"]/summaryStats["VUS_syn","sd"],
            `Antagonism Score` = summaryStats["VUS_ant","mean"]/summaryStats["VUS_ant","sd"],
            # Finally some information about model fit
            `s` = summaryStats["s", "mean"],
            `returnCode` = fit$returnCode,
            `divergentTransitions` = sum(fit$divergent),

            check.names = FALSE, stringsAsFactors = FALSE)
          if (fit$model$bayes_factor){
            synMetrics$`Bayes Factor` = fit$bayesfactor
          }

        } else {
          # Create some summaries
          summaryStats = rstan::summary(fit$stanfit)$summary
          synMetrics <- data.frame(
            `listID` = i,
            `Experiment ID` = data$experiment_ID,
            `Drug A` = data$drug_names[1],
            `Drug B` = data$drug_names[2],
            # Include some summary stats for the experiment
            `EC50 (Drug A)` = NA,
            `EC50 (Drug B)` = NA,
            `DSS (Drug A)` = NA,
            `DSS (Drug B)` = NA,
            `Synergy (mean)` = NA,
            `Synergy (sd)` = NA,
            `Antagonism (mean)` = NA,
            `Antagonism (sd)` = NA,
            `Efficacy (mean)` = NA,
            `Efficacy (sd)` = NA,
            `Non-interaction Efficacy (mean)` = NA,
            `Non-interaction Efficacy (sd)` = NA,
            `Interaction (mean)` = NA,
            `Interaction (sd)` = NA,

            # Calculate standardized scores for comparisons
            `Synergy Score` = NA,
            `Antagonism Score` = NA,
            # Finally some information about model fit
            `s` = NA,
            `returnCode` = 2,
            `divergentTransitions` = NA,

            check.names = FALSE, stringsAsFactors = FALSE)
          if (fit$model$bayes_factor){
            synMetrics$`Bayes Factor` = NA
          }
        }
      })

      # Garbage collection
      dir <- tempdir()
      file.names <- list.files(path = dir, pattern="*.csv")
      unlink(paste0(dir,"/",file.names))

      #Return
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
  # Set those with too high noise as returnCode 2
  toReturn$screenSummary[which(toReturn$screenSummary[,"s"] > 1), "returnCode"] = 2

  # Then remove all with returnCode 2
  toReturn$screenSummary = toReturn$screenSummary[which(!(toReturn$screenSummary[,"returnCode"] == 2)),]

  # Do some checks to see if things converges okay

  propFailed = mean((toReturn$screenSummary[,"returnCode"]>0),na.rm = T)

  if (propFailed > 0){
    if (bayesynergy_params$method == "sampling"){
      warning(paste0(round(100*propFailed,digits=2),"% of experiments returned a warning or error -- check these."))
    }
    if (bayesynergy_params$method == "vb"){
      warning(paste0(round(100*propFailed,digits=2),"% of experiments returned a warning or error. NOTE: For method = 'vb' this is to be expected"))
    }
  }
  if (sum((toReturn$screenSummary[,"returnCode"]!=2)) != length(experiments)){
    diff = length(experiments) - sum((toReturn$screenSummary[,"returnCode"]!=2))
    warning(paste(diff," experiment(s) failed to process"))
  }

  # Pick out failed experiments and put in list
  failedID = unique(c(setdiff(1:length(experiments),toReturn$screenSummary$`listID`),which(toReturn$screenSummary[,"returnCode"]==2)))
  failedExperiments = experiments[failedID]
  if (length(failedID>0)){
    toReturn$failed = failedExperiments
  }

  # Finally, if we are re-running, the new runs are combined with the old
  if (rerun){
    toReturn$screenSummary = rbind(oldrun,toReturn$screenSummary)
  }

  # Returning
  class(toReturn) <- "synergyscreen"
  return(toReturn)

}

