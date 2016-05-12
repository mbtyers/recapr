# recapr 

### Estimating, Testing, and Simulating Abundance in a Mark-Recapture

Tools are provided for estimating, testing, and simulating abundance in a two-event (Petersen) mark-recapture experiment.  Functions are given to calculate the Petersen, Chapman, and Bailey estimators and associated variances.  However, the principal utility is a set of functions to simulate random draws from these estimators, and use these to conduct hypothesis tests and power calculations.  Additionally, a set of functions are provided for generating confidence intervals via bootstrapping.

### Commonly-used functions


### Installation

The 'riverdist' package is currently available on Github, and can be installed in R with the following code:

`install.packages("devtools",dependencies=T)`

`devtools::install_github("mbtyers/recapr")`

### Issues

This package has no known issues.  

[![Travis-CI Build Status](https://travis-ci.org/mbtyers/recapr.svg?branch=master)](https://travis-ci.org/mbtyers/recapr)
