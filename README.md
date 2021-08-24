# bayesynergy: flexible Bayesian modelling of synergistic interaction effects in in-vitro drug combination experiments

An R package for Bayesian semi-parametric modelling of in-vitro drug combination experiments 

[![DOI](https://zenodo.org/badge/192396008.svg)](https://zenodo.org/badge/latestdoi/192396008)
[![Build Status](https://travis-ci.org/ocbe-uio/bayesynergy.svg?branch=master)](https://travis-ci.org/ocbe-uio/bayesynergy)



<img src="https://github.com/ocbe-uio/bayesynergy/blob/master/doc/workflow.png?raw=true" width="90%" align="center"></img>

The bayesynergy package implements a Bayesian semi-parametric model for the drug combination experiment. Through an efficient implementation in [Stan](https://mc-stan.org/users/interfaces/rstan), the model provides estimates of the full dose-response surface with uncertainty quantification. From the posterior dose-response function, estimates of synergy and antagonism are derived that better reflect the true uncertainty in these estimates. The model handles incomplete and messy datasets, naturally includes replicates, and contains parallel processing for large drug combination screens.

[Paper](https://doi.org/10.1093/bib/bbab251)

## Installation
Prior to installing the package, please install [RStan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started), and make sure you can run the models there. Stan requires the configuration of a C++ toolchain, which can be tricky on some systems. 


```r
install.packages('devtools')
library(devtools)
install_github('ocbe-uio/bayesynergy', build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
```

## Citation


Leiv Rønneberg, Andrea Cremaschi, Robert Hanes, Jorrit M Enserink, Manuela Zucknick, bayesynergy: flexible Bayesian modelling of synergistic interaction effects in <em>in vitro</em> drug combination experiments, <em>Briefings in Bioinformatics</em>, 2021;, bbab251, https://doi.org/10.1093/bib/bbab251



## Updates
- 26-05-21: Updated to version 2.4.1, minor notation changes. Ability to calculate the Bayes Factor comparing models with and without interaction
- 29-03-21: Updated to version 2.4, some minor tweaks to priors, default settings and plotting function. Added option to set hyperparameters of observation noise
- 13-01-21: Updated to version 2.3, added new dataset for use with synergyscreen(), fixed bug with synergyscreen() filling up hard drive, added vignette
- 01-12-20: Updated to version 2.2, added additional plotting functions for the synergyscreen() function, added a heteroscedastic setting for the noise term, + some minor improvements overall.
- 14-10-20: Updated to version 2.1, removed Gibbs sampler, should now compile easier on Windows systems
- 27-08-20: Updated to version 2.0, now running on Stan!
- 03-03-20: Added support for missing values!
- 14-04-20: Added the option to estimate lower asymptotes for monotherapies
