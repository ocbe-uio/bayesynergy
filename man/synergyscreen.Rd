% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synergyscreen.R
\name{synergyscreen}
\alias{synergyscreen}
\title{Function for fitting high-througput drug combination screens with parallel processing}
\usage{
synergyscreen(
  experiments,
  return_samples = F,
  save_raw = T,
  save_plots = T,
  path = NULL,
  parallel = T,
  max_cores = NULL,
  max_retries = 3,
  plot_params = list(),
  bayesynergy_params = list()
)
}
\arguments{
\item{experiments}{A list of experiments obtained from a high-throughput screen. See *Details* for more information on the structure of each element.}

\item{return_samples}{logical; if TRUE, the function returns the full fitted \code{\link[bayesynergy]{bayesynergy}} object.}

\item{save_raw}{logical; if TRUE, the raw bayesynergy object is saved for each individual experiment.}

\item{save_plots}{logical; if TRUE, plots for each individual experiment is saved.}

\item{path}{string; path for saving output and plot for each individual experiment.}

\item{parallel}{logical; if TRUE, parallel processing is utilized to run the screen.}

\item{max_cores}{integer; the maximum number of cores to utilize for the parallel processing.}

\item{max_retries}{integer; the maximum number of retries utilized in model fit.}

\item{plot_params}{list; parameters to be passed to the plotting function. See \link{plot.bayesynergy} for details.}

\item{bayesynergy_params}{list; parameters to be passed to the bayesynergy function. See \link{bayesynergy} for details.}
}
\value{
A list containing two elements
\tabular{ll}{
screenSummary \tab data frame, posterior summary statistics for each experiment. \cr
failed \tab A list containing elements experiments that failed to process. \cr
screenSamples \tab if requested, a list containing the fitted \code{\link[bayesynergy]{bayesynergy}} object for each experiment.
}
}
\description{
The function \code{synergyscreen} allows the fitting of high-throughput drug combination screens through parallel processing
}
\details{
The elements of \strong{experiments} must themselves be lists with the following elements
\tabular{ll}{
y: \tab matrix of viability scores \cr
x: \tab matrix of concentrations \cr
drug_names: \tab vector of drug names in the experiment \cr
experiment_ID: \tab string denoting the unique experiment ID, e.g. cell line name.
units vector of size 2; concentration units for the drugs, e.g. c("\eqn{\mu}M","\eqn{\mu}M")
}
}
\examples{
\dontrun{
library(bayesynergy)
data("mathews_DLBCL")
experiment1 = list(y = mathews_DLBCL[[1]][[1]],
x = mathews_DLBCL[[1]][[2]],
drug_names = c("ispinesib","ibrutinib"))
experiment2 = list(y = mathews_DLBCL[[2]][[1]],
x = mathews_DLBCL[[2]][[2]],
drug_names = c("canertinib","ibrutinib"))
experiments = list(experiment1,experiment2)
fit <- synergyscreen(experiments)
}

}
