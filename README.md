# BayeSyneRgy

An R package for Bayesian semi-parametric modelling of in-vitro drug combination experiments 

[![DOI](https://zenodo.org/badge/192396008.svg)](https://zenodo.org/badge/latestdoi/192396008)
[![Build Status](https://travis-ci.org/ocbe-uio/bayesynergy.svg?branch=master)](https://travis-ci.org/ocbe-uio/bayesynergy)

This package rewrites and extends the functionalities of the [synergysplines](https://github.com/ocbe-uio/synergysplines) package, offering significant computational gains and other fixes. More information [here](https://github.com/ocbe-uio/bayesynergy/issues/2#issuecomment-718553824).

## Updates
- 14-10-20: Updated to version 2.1, removed Gibbs sampler, should now compile easier on Windows systems
- 27-08-20: Updated to version 2.0, now running on Stan!
- 03-03-20: Added support for missing values!
- 14-04-20: Added the option to estimate lower asymptotes for monotherapies



## Installation

```r
install.packages('devtools')
library(devtools)
install_github('ocbe-uio/bayesynergy',build_opts = c("--no-resave-data", "--no-manual"))
```
