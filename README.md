# recapr 

### Estimating, Testing, and Simulating Abundance in a Mark-Recapture

Tools are provided for estimating, testing, and simulating abundance in a two-event (Petersen) mark-recapture experiment.  Functions are given to calculate the Petersen, Chapman, and Bailey estimators and associated variances.  However, the principal utility is a set of functions to simulate random draws from these estimators, and use these to conduct hypothesis tests and power calculations.  Additionally, a set of functions are provided for generating confidence intervals via bootstrapping.

### Commonly-used functions

* `NChapman()`, `NPetersen()`, and `NBailey()` calculate the values of Chapman, Petersen, or Bailey abundance estimates, given values of sample sizes and number of recaptures 

* `vChapman()`, `vPetersen()`, and `vBailey()` calculate the estimated variance of Chapman, Petersen, or Bailey abundance estimates, given values of sample sizes and number of recaptures, and `seChapman()`, `sePetersen()`, and `seBailey()` give standard errors

* `rChapman()`, `rPetersen()`, and `rBailey()` return vectors of random draws from the Chapman, Petersen, or Bailey abundance estimates, given values of true abundance and sample sizes

* `pChapman()`, `pPetersen()`, and `pBailey()` use many random draws to calculate approximate p-values for hypothesis testing

* `powChapman()`, `powPetersen()`, and `powBailey()` use simulation to calculate hypothesis testing power, given alternative abundance

* `ciChapman()`, `ciPetersen()`, and `ciBailey()` calculate confidence intervals for abundance using bootstrapping and/or normal approximation

* `plotdiscdensity()` produces an empirical pmf plot of a vector of discrete values, such as that returned from an abundance estimate simulation, that is more appropriate than a traditional kernel density plot and perhaps more illustrative than a histogram

### Installation

The 'recapr' package is currently available on Github, and can be installed in R with the following code:

`install.packages("devtools",dependencies=T)`

`devtools::install_github("mbtyers/recapr")`

### Issues

This package has no known issues.  

[![Travis-CI Build Status](https://travis-ci.org/mbtyers/recapr.svg?branch=master)](https://travis-ci.org/mbtyers/recapr)
