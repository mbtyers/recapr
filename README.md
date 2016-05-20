# recapr 

### Estimating, Testing, and Simulating Abundance in a Mark-Recapture

Tools are provided for estimating, testing, and simulating abundance in a two-event (Petersen) mark-recapture experiment.  Functions are given to calculate the Petersen, Chapman, and Bailey estimators and associated variances.  However, the principal utility is a set of functions to simulate random draws from these estimators, and use these to conduct hypothesis tests and power calculations.  Additionally, a set of functions are provided for generating confidence intervals via bootstrapping.  Functions are also provided to test abundance estimator consistency under complete or partial stratification, and to calculate stratified or Darroch estimators.  Functions are also provided to calculate recommended sample sizes.

### Commonly-used functions

* `NChapman()`, `NPetersen()`, and `NBailey()` calculate the values of Chapman, Petersen, or Bailey abundance estimates, given values of sample sizes and number of recaptures 

* `vChapman()`, `vPetersen()`, and `vBailey()` calculate the estimated variance of Chapman, Petersen, or Bailey abundance estimates, given values of sample sizes and number of recaptures, and `seChapman()`, `sePetersen()`, and `seBailey()` give standard errors

* `rChapman()`, `rPetersen()`, and `rBailey()` return vectors of random draws from the Chapman, Petersen, or Bailey abundance estimates, given values of true abundance and sample sizes

* `pChapman()`, `pPetersen()`, and `pBailey()` use many random draws to calculate approximate p-values for hypothesis testing

* `powChapman()`, `powPetersen()`, and `powBailey()` use simulation to calculate hypothesis testing power, given alternative abundance

* `ciChapman()`, `ciPetersen()`, and `ciBailey()` calculate confidence intervals for abundance using bootstrapping and/or normal approximation

* `plotdiscdensity()` produces an empirical pmf plot of a vector of discrete values, such as that returned from an abundance estimate simulation, that is more appropriate than a traditional kernel density plot and perhaps more illustrative than a histogram

* `consistencytest()` and `strattest()` provide the typical chi-squared tests for the consistency of a Petersen-type estimator, and provide evidence of the necessity of a stratified or partially stratified (Darroch-type) estimator

* `Nstrat()`, `vstrat()`, `sestrat()` and `cistrat()` provide estimation if a completely stratified estimator is used

* `NDarroch()` provides estimation if a spatially or temporally stratified estimator is used, or if strata differs between sampling events

* `n2RR()` provides recommended sample size using Robson-Regier, and `plotn2sim()` and `plotn1n2simmatrix()` provide graphical explorations of recommended sample sizes via simulation

### Installation

The 'recapr' package is currently available on Github, and can be installed in R with the following code:

`install.packages("devtools",dependencies=T)`

`devtools::install_github("mbtyers/recapr")`

### Issues

This package has no known issues.  

[![Travis-CI Build Status](https://travis-ci.org/mbtyers/recapr.svg?branch=master)](https://travis-ci.org/mbtyers/recapr)
