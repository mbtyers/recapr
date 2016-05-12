#' Confidence Intervals for the Petersen Estimator
#' @description Calculates approximate confidence intervals(s) for the Petersen
#'   estimator, using bootstrapping, the Normal approximation, or both.
#'
#'   The bootstrap interval is created by resampling the data in the second
#'   sampling event, with replacement; that is, drawing bootstrap values of m2
#'   from a binomial distribution with probability parameter m2/n2.  This
#'   technique has been shown to better approximate the distribution of the
#'   abundance estimator.  Resulting CI endpoints both have larger values than
#'   those calculated from a normal distribution, but this better captures the
#'   positive skew of the estimator.  Coverage has been investigated by means of
#'   simulation under numerous scenarios and has consistently outperformed the
#'   normal interval.  The user is welcomed to investigate the coverage under
#'   relevant scenarios.
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of marked individuals recaptured in the second sample
#' @param conf The confidence level of the desired intervals.  Defaults to 0.95.
#' @param method Which method of confidence interval to return.  Allowed values
#'   are \code{"norm"}, \code{"boot"}, or \code{"both"}.  Defaults to
#'   \code{"both"}.
#' @param bootreps Number of bootstrap replicates to use.  Defaults to 10000.
#' @param useChapvar Whether to use the Chapman estimator variance instead of
#'   the Petersen estimator variance for the normal-distribution interval.
#'   Defaults to \code{FALSE}.
#' @return A list with the abundance estimate and confidence interval bounds for
#'   the normal-distribution and/or bootstrap confidence intervals.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NPetersen}, \link{vPetersen}, \link{sePetersen},
#'   \link{rPetersen}, \link{pPetersen}, \link{powPetersen}
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats quantile
#' @examples
#' ciPetersen(n1=100, n2=100, m2=20)
#' @export
ciPetersen <- function(n1, n2, m2, conf=0.95, method="both", bootreps=10000,useChapvar=FALSE) {
  if(!any(method==c("boot","norm","both"))) stop("invalid method specified")
  bounds <- c((1-conf)/2,1-((1-conf)/2))
  out <- list()
  out$Nhat <- NPetersen(n1, n2, m2)
  if(method=="norm" | method=="both") {
    se <- ifelse(useChapvar,seChapman(n1, n2, m2), sePetersen(n1, n2, m2))
    ciNorm <- out$Nhat + qnorm(bounds)*se
    out$ciNorm <- ciNorm
  }
  if(method=="boot" | method=="both") {
    m2boot <- rbinom(bootreps,size=n2,prob=(m2/n2))
    Nhatboot <- NPetersen(n1,n2,m2boot)
    out$ciBoot <- unname(quantile(Nhatboot,bounds,na.rm=T))
  }
  return(out)
}


#' Confidence Intervals for the Chapman Estimator
#' @description Calculates approximate confidence intervals(s) for the Chapman
#'   estimator, using bootstrapping, the Normal approximation, or both.
#'
#'   The bootstrap interval is created by resampling the data in the second
#'   sampling event, with replacement; that is, drawing bootstrap values of m2
#'   from a binomial distribution with probability parameter m2/n2.  This
#'   technique has been shown to better approximate the distribution of the
#'   abundance estimator.  Resulting CI endpoints both have larger values than
#'   those calculated from a normal distribution, but this better captures the
#'   positive skew of the estimator.  Coverage has been investigated by means of
#'   simulation under numerous scenarios and has consistently outperformed the
#'   normal interval.  The user is welcomed to investigate the coverage under
#'   relevant scenarios.
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of marked individuals recaptured in the second sample
#' @param conf The confidence level of the desired intervals.  Defaults to 0.95.
#' @param method Which method of confidence interval to return.  Allowed values
#'   are \code{"norm"}, \code{"boot"}, or \code{"both"}.  Defaults to
#'   \code{"both"}.
#' @param bootreps Number of bootstrap replicates to use.  Defaults to 10000.
#' @return A list with the abundance estimate and confidence interval bounds for
#'   the normal-distribution and/or bootstrap confidence intervals.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{vChapman}, \link{seChapman},
#'   \link{rChapman}, \link{pChapman}, \link{powChapman}
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats quantile
#' @examples
#' ciChapman(n1=100, n2=100, m2=20)
#' @export
ciChapman <- function(n1, n2, m2, conf=0.95, method="both", bootreps=10000) {
  if(!any(method==c("boot","norm","both"))) stop("invalid method specified")
  bounds <- c((1-conf)/2,1-((1-conf)/2))
  out <- list()
  out$Nhat <- NChapman(n1, n2, m2)
  if(method=="norm" | method=="both") {
    se <- seChapman(n1, n2, m2)
    ciNorm <- out$Nhat + qnorm(bounds)*se
    out$ciNorm <- ciNorm
  }
  if(method=="boot" | method=="both") {
    m2boot <- rbinom(bootreps,size=n2,prob=(m2/n2))
    Nhatboot <- NChapman(n1,n2,m2boot)
    out$ciBoot <- unname(quantile(Nhatboot,bounds,na.rm=T))
  }
  return(out)
}


#' Confidence Intervals for the Bailey Estimator
#' @description Calculates approximate confidence intervals(s) for the Bailey
#'   estimator, using bootstrapping, the Normal approximation, or both.
#'
#'   The bootstrap interval is created by resampling the data in the second
#'   sampling event, with replacement; that is, drawing bootstrap values of m2
#'   from a binomial distribution with probability parameter m2/n2.  This
#'   technique has been shown to better approximate the distribution of the
#'   abundance estimator.  Resulting CI endpoints both have larger values than
#'   those calculated from a normal distribution, but this better captures the
#'   positive skew of the estimator.  Coverage has been investigated by means of
#'   simulation under numerous scenarios and has consistently outperformed the
#'   normal interval.  The user is welcomed to investigate the coverage under
#'   relevant scenarios.
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of marked individuals recaptured in the second sample
#' @param conf The confidence level of the desired intervals.  Defaults to 0.95.
#' @param method Which method of confidence interval to return.  Allowed values
#'   are \code{"norm"}, \code{"boot"}, or \code{"both"}.  Defaults to
#'   \code{"both"}.
#' @param bootreps Number of bootstrap replicates to use.  Defaults to 10000.
#' @return A list with the abundance estimate and confidence interval bounds for
#'   the normal-distribution and/or bootstrap confidence intervals.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NBailey}, \link{vBailey}, \link{seBailey},
#'   \link{rBailey}, \link{pBailey}, \link{powBailey}
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats quantile
#' @examples
#' ciBailey(n1=100, n2=100, m2=20)
#' @export
ciBailey <- function(n1, n2, m2, conf=0.95, method="both", bootreps=10000) {
  if(!any(method==c("boot","norm","both"))) stop("invalid method specified")
  bounds <- c((1-conf)/2,1-((1-conf)/2))
  out <- list()
  out$Nhat <- NBailey(n1, n2, m2)
  if(method=="norm" | method=="both") {
    se <- seBailey(n1, n2, m2)
    ciNorm <- out$Nhat + qnorm(bounds)*se
    out$ciNorm <- ciNorm
  }
  if(method=="boot" | method=="both") {
    m2boot <- rbinom(bootreps,size=n2,prob=(m2/n2))
    Nhatboot <- NBailey(n1,n2,m2boot)
    out$ciBoot <- unname(quantile(Nhatboot,bounds,na.rm=T))
  }
  return(out)
}

#ciPetersen(100,100,10)
#
# nsim <- 10000
# n1 <- 300
# n2 <- 200
# m2 <- 50
# N <- 1200
# m2draws <- rhyper(nsim,n1,N-n1,n2)
# covnorm <- covboot <- NA
# for(i in 1:nsim) {
#   cis <- ciPetersen(n1,n2,m2draws[i],conf=.9,useChapvar = T)
#   covnorm[i] <- (cis$ciNorm[1] <= N) & (cis$ciNorm[2] >= N)
#   covboot[i] <- (cis$ciBoot[1] <= N) & (cis$ciBoot[2] >= N)
# }
# mean(covnorm)
# mean(covboot)
#
# mean((rChapman(10000,N,n1,n2)))
