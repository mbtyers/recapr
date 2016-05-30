#' Random Draws from the Chapman Estimator
#' @description Returns a vector of random draws from the Chapman estimator in a
#'   mark-recapture experiment, given values of the true abundance and the
#'   sample size in both events.  The function first simulates a vector of
#'   recaptures (m2) from a hypergeometric distribution, and then uses these to
#'   compute a vector of draws from the estimator.
#'
#' If capture probabilities (\code{p1} and/or \code{p2}) are specified instead of sample size(s), the sample size(s) will first be drawn from a binomial distribution, then the number of recaptures.  If both sample size and capture probability are specified for a given sampling event, only the sample size will be used.
#' @param length The length of the random vector to return.
#' @param N The value of the true abundance. This may be a single number or
#'   vector of values equal to \code{length}.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values equal to \code{length}.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values equal to \code{length}.
#' @param p1 Alternately, probability of capture in the first sample. This
#'   may be a single number or vector of values equal to \code{length}.
#' @param p2 Alternately, probability of capture in the second sample.  This may be a
#'   single number or vector of values equal to \code{length}.
#' @return A vector of random draws from the Chapman estimator
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{vChapman}, \link{seChapman}, \link{pChapman},
#'   \link{powChapman}, \link{ciChapman}
#' @importFrom stats rhyper
#' @importFrom stats rbinom
#' @examples
#' draws <- rChapman(length=100000, N=500, n1=100, n2=100)
#' plotdiscdensity(draws)  #plots the density of a vector of discrete values
#' @export
rChapman <- function(length,N,n1=NULL,n2=NULL,p1=NULL,p2=NULL) {
  if(length(N)!=1 & length(N!=length)) stop("N must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n1)) if(length(n1)!=1 & length(n1)!=length) stop("n1 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n2)) if(length(n2)!=1 & length(n2)!=length) stop("n2 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(p1)) if(length(p1)!=1 & length(p1)!=length) stop("p1 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(p2)) if(length(p2)!=1 & length(p2)!=length) stop("p2 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n1) & !is.null(p1)) warning("both n1 and p1 specified - only n1 used")
  if(!is.null(n2) & !is.null(p2)) warning("both n2 and p2 specified - only n2 used")
  if(is.null(n1)) {
    if(is.null(p1)) stop("need to supply either n1 or p1")
    n1 <- rbinom(n=length, size=N, prob=p1)    # this is new
  }
  if(is.null(n2)) {
    if(is.null(p2)) stop("need to supply either n2 or p2")
    n2 <- rbinom(n=length, size=N, prob=p2)    # this is new
  }
  m2 <- rhyper(length, n1, N-n1, n2)
  Nhat <- (n1+1)*(n2+1)/(m2+1)-1
  return(Nhat)
}


#' Random Draws from the Petersen Estimator
#' @description Returns a vector of random draws from the Petersen estimator in a
#'   mark-recapture experiment, given values of the true abundance and the
#'   sample size in both events.  The function first simulates a vector of
#'   recaptures (m2) from a hypergeometric distribution, and then uses these to
#'   compute a vector of draws from the estimator.
#'
#' If capture probabilities (\code{p1} and/or \code{p2}) are specified instead of sample size(s), the sample size(s) will first be drawn from a binomial distribution, then the number of recaptures.  If both sample size and capture probability are specified for a given sampling event, only the sample size will be used.
#' @param length The length of the random vector to return.
#' @param N The value of the true abundance. This may be a single number or
#'   vector of values equal to \code{length}.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values equal to \code{length}.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values equal to \code{length}.
#' @param p1 Alternately, probability of capture in the first sample. This
#'   may be a single number or vector of values equal to \code{length}.
#' @param p2 Alternately, probability of capture in the second sample.  This may be a
#'   single number or vector of values equal to \code{length}.
#' @return A vector of random draws from the Petersen estimator
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NPetersen}, \link{vPetersen}, \link{sePetersen}, \link{pPetersen},
#'   \link{powPetersen}, \link{ciPetersen}
#' @importFrom stats rhyper
#' @importFrom stats rbinom
#' @examples
#' draws <- rPetersen(length=100000, N=500, n1=100, n2=100)
#' plotdiscdensity(draws)  #plots the density of a vector of discrete values
#' @export
rPetersen <- function(length,N,n1=NULL,n2=NULL,p1=NULL,p2=NULL) {
  if(length(N)!=1 & length(N!=length)) stop("N must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n1)) if(length(n1)!=1 & length(n1)!=length) stop("n1 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n2)) if(length(n2)!=1 & length(n2)!=length) stop("n2 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(p1)) if(length(p1)!=1 & length(p1)!=length) stop("p1 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(p2)) if(length(p2)!=1 & length(p2)!=length) stop("p2 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n1) & !is.null(p1)) warning("both n1 and p1 specified - only n1 used")
  if(!is.null(n2) & !is.null(p2)) warning("both n2 and p2 specified - only n2 used")
  if(is.null(n1)) {
    if(is.null(p1)) stop("need to supply either n1 or p1")
    n1 <- rbinom(n=length, size=N, prob=p1)    # this is new
  }
  if(is.null(n2)) {
    if(is.null(p2)) stop("need to supply either n2 or p2")
    n2 <- rbinom(n=length, size=N, prob=p2)    # this is new
  }
  m2 <- rhyper(length, n1, N-n1, n2)
  Nhat <- (n1)*(n2)/(m2)
  return(Nhat)
}


#' Random Draws from the Bailey Estimator
#' @description Returns a vector of random draws from the Bailey estimator in a
#'   mark-recapture experiment, given values of the true abundance and the
#'   sample size in both events.  The function first simulates a vector of
#'   recaptures (m2) from a binomial distribution, and then uses these to
#'   compute a vector of draws from the estimator.
#'
#' If capture probabilities (\code{p1} and/or \code{p2}) are specified instead of sample size(s), the sample size(s) will first be drawn from a binomial distribution, then the number of recaptures.  If both sample size and capture probability are specified for a given sampling event, only the sample size will be used.
#' @param length The length of the random vector to return.
#' @param N The value of the true abundance. This may be a single number or
#'   vector of values equal to \code{length}.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values equal to \code{length}.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values equal to \code{length}.
#' @param p1 Alternately, probability of capture in the first sample. This
#'   may be a single number or vector of values equal to \code{length}.
#' @param p2 Alternately, probability of capture in the second sample.  This may be a
#'   single number or vector of values equal to \code{length}.
#' @return A vector of random draws from the Bailey estimator
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NBailey}, \link{vBailey}, \link{seBailey}, \link{pBailey},
#'   \link{powBailey}, \link{ciBailey}
#' @importFrom stats rbinom
#' @examples
#' draws <- rBailey(length=100000, N=500, n1=100, n2=100)
#' plotdiscdensity(draws)  #plots the density of a vector of discrete values
#' @export
rBailey <- function(length,N,n1=NULL,n2=NULL,p1=NULL,p2=NULL) {
  if(length(N)!=1 & length(N!=length)) stop("N must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n1)) if(length(n1)!=1 & length(n1)!=length) stop("n1 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n2)) if(length(n2)!=1 & length(n2)!=length) stop("n2 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(p1)) if(length(p1)!=1 & length(p1)!=length) stop("p1 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(p2)) if(length(p2)!=1 & length(p2)!=length) stop("p2 must be a single number, or a vector of length equal to the number of draws")
  if(!is.null(n1) & !is.null(p1)) warning("both n1 and p1 specified - only n1 used")
  if(!is.null(n2) & !is.null(p2)) warning("both n2 and p2 specified - only n2 used")
  if(is.null(n1)) {
    if(is.null(p1)) stop("need to supply either n1 or p1")
    n1 <- rbinom(n=length, size=N, prob=p1)    # this is new
  }
  if(is.null(n2)) {
    if(is.null(p2)) stop("need to supply either n2 or p2")
    n2 <- rbinom(n=length, size=N, prob=p2)    # this is new
  }
  #m2 <- rhyper(length, n1, N-n1, n2)
  m2 <- rbinom(length,n2,n1/N)
  Nhat <- (n1)*(n2+1)/(m2+1)
  return(Nhat)
}


#' Hypothesis Testing Using the Bailey Estimator
#' @description Approximates a p-value for a hypothesis test of the Bailey
#'   estimator by means of many simulated draws from the null distribution, conditioned on sample sizes.
#' @param estN The estimated abundance.  Either this or the number of recaptures
#'   (\code{m2}) must be specified.
#' @param nullN The abundance given by the null hypothesis
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of recaptures.  Either this or the estimated abundance
#'   (\code{estN}) must be specified.
#' @param nsim Number of simulated values to draw.  Defaults to 100000.
#' @param alternative Direction of the alternative hypothesis.  Allowed values
#'   are \code{"less"}, \code{"greater"}, or \code{"2-sided"}.  Defaults to \code{"less"}.
#' @return An approximate p-value for the specified hypothesis test.  If
#'   \code{m2} is specified rather than \code{estN}, output will be returned as
#'   a list with two elements: the estimated abundance and p-value.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NBailey}, \link{vBailey}, \link{seBailey}, \link{rBailey},
#'   \link{powBailey}, \link{ciBailey}
#' @examples
#' output <- pBailey(nullN=500, n1=100, n2=100, m2=28)
#' output
#'
#' plotdiscdensity(rBailey(length=100000, N=500, n1=100, n2=100))
#' abline(v=output$estN, lwd=2, col=2)
#' abline(v=500, lwd=2, lty=2)
#'
#'
#' output <- pBailey(nullN=500, n1=100, n2=100, m2=28, alternative="2-sided")
#' output
#'
#' plotdiscdensity(rBailey(length=100000, N=500, n1=100, n2=100))
#' twosided <- 500 + c(-1,1)*abs(500-output$estN)
#' abline(v=twosided, lwd=2, col=2)
#' abline(v=500, lwd=2, lty=2)
#' @export
pBailey <- function(estN=NULL,nullN,n1,n2,m2=NULL,nsim=100000,alternative="less") {
  if(!any(alternative==c("less","greater","2-sided"))) stop("invalid value for alternative hypothesis")
  printest <- F
  if(is.null(estN) & !is.null(m2)) {
    estN <- NBailey(n1=n1,n2=n2,m2=m2)
    printest <- T
  }
  if(is.null(estN) & is.null(m2)) stop("estN or m2 must be provided")
  asdf <- rBailey(length=nsim,N=nullN,n1=n1,n2=n2)
  if(alternative=="greater") pval <- mean(asdf>=estN)
  if(alternative=="less") pval <- mean(asdf<=estN)
  if(alternative=="2-sided") pval <- mean(abs(asdf-nullN)>=abs(estN-nullN))        # this may not be right
  if(printest) {
    out <- list(estN=estN,pval=pval)
  }
  if(!printest) out <- pval
  return(out)
}


#' Hypothesis Testing Using the Chapman Estimator
#' @description Approximates a p-value for a hypothesis test of the Chapman
#'   estimator by means of many simulated draws from the null distribution, conditioned on sample sizes.
#' @param estN The estimated abundance.  Either this or the number of recaptures
#'   (\code{m2}) must be specified.
#' @param nullN The abundance given by the null hypothesis
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of recaptures.  Either this or the estimated abundance
#'   (\code{estN}) must be specified.
#' @param nsim Number of simulated values to draw.  Defaults to 100000.
#' @param alternative Direction of the alternative hypothesis.  Allowed values
#'   are \code{"less"}, \code{"greater"}, or \code{"2-sided"}.  Defaults to \code{"less"}.
#' @return An approximate p-value for the specified hypothesis test.  If
#'   \code{m2} is specified rather than \code{estN}, output will be returned as
#'   a list with two elements: the estimated abundance and p-value.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{vChapman}, \link{seChapman}, \link{rChapman},
#'   \link{powChapman}, \link{ciChapman}
#' @examples
#' output <- pChapman(nullN=500, n1=100, n2=100, m2=28)
#' output
#'
#' plotdiscdensity(rChapman(length=100000, N=500, n1=100, n2=100))
#' abline(v=output$estN, lwd=2, col=2)
#' abline(v=500, lwd=2, lty=2)
#'
#'
#' output <- pChapman(nullN=500, n1=100, n2=100, m2=28, alternative="2-sided")
#' output
#'
#' plotdiscdensity(rChapman(length=100000, N=500, n1=100, n2=100))
#' twosided <- 500 + c(-1,1)*abs(500-output$estN)
#' abline(v=twosided, lwd=2, col=2)
#' abline(v=500, lwd=2, lty=2)
#' @export
pChapman <- function(estN=NULL,nullN,n1,n2,m2=NULL,nsim=100000,alternative="less") {
  if(!any(alternative==c("less","greater","2-sided"))) stop("invalid value for alternative hypothesis")
  printest <- F
  if(is.null(estN) & !is.null(m2)) {
    estN <- NChapman(n1=n1,n2=n2,m2=m2)
    printest <- T
  }
  if(is.null(estN) & is.null(m2)) stop("estN or m2 must be provided")
  asdf <- rChapman(length=nsim,N=nullN,n1=n1,n2=n2)
  if(alternative=="greater") pval <- mean(asdf>=estN)
  if(alternative=="less") pval <- mean(asdf<=estN)
  if(alternative=="2-sided") pval <- mean(abs(asdf-nullN)>=abs(estN-nullN))        # this may not be right
  if(printest) {
    out <- list(estN=estN,pval=pval)
  }
  if(!printest) out <- pval
  return(out)
}


#' Hypothesis Testing Using the Petersen Estimator
#' @description Approximates a p-value for a hypothesis test of the Petersen
#'   estimator by means of many simulated draws from the null distribution, conditioned on sample sizes.
#' @param estN The estimated abundance.  Either this or the number of recaptures
#'   (\code{m2}) must be specified.
#' @param nullN The abundance given by the null hypothesis
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of recaptures.  Either this or the estimated abundance
#'   (\code{estN}) must be specified.
#' @param nsim Number of simulated values to draw.  Defaults to 100000.
#' @param alternative Direction of the alternative hypothesis.  Allowed values
#'   are \code{"less"}, \code{"greater"}, or \code{"2-sided"}.  Defaults to \code{"less"}.
#' @return An approximate p-value for the specified hypothesis test.  If
#'   \code{m2} is specified rather than \code{estN}, output will be returned as
#'   a list with two elements: the estimated abundance and p-value.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NPetersen}, \link{vPetersen}, \link{sePetersen}, \link{rPetersen},
#'   \link{powPetersen}, \link{ciPetersen}
#' @examples
#' output <- pPetersen(nullN=500, n1=100, n2=100, m2=28)
#' output
#'
#' plotdiscdensity(rPetersen(length=100000, N=500, n1=100, n2=100))
#' abline(v=output$estN, lwd=2, col=2)
#' abline(v=500, lwd=2, lty=2)
#'
#'
#' output <- pPetersen(nullN=500, n1=100, n2=100, m2=28, alternative="2-sided")
#' output
#'
#' plotdiscdensity(rPetersen(length=100000, N=500, n1=100, n2=100))
#' twosided <- 500 + c(-1,1)*abs(500-output$estN)
#' abline(v=twosided, lwd=2, col=2)
#' abline(v=500, lwd=2, lty=2)
#' @export
pPetersen <- function(estN=NULL,nullN,n1,n2,m2=NULL,nsim=100000,alternative="less") {
  if(!any(alternative==c("less","greater","2-sided"))) stop("invalid value for alternative hypothesis")
  printest <- F
  if(is.null(estN) & !is.null(m2)) {
    estN <- NPetersen(n1=n1,n2=n2,m2=m2)
    printest <- T
  }
  if(is.null(estN) & is.null(m2)) stop("estN or m2 must be provided")
  asdf <- rPetersen(length=nsim,N=nullN,n1=n1,n2=n2)
  if(alternative=="greater") pval <- mean(asdf>=estN)
  if(alternative=="less") pval <- mean(asdf<=estN)
  if(alternative=="2-sided") pval <- mean(abs(asdf-nullN)>=abs(estN-nullN))        # this may not be right
  if(printest) {
    out <- list(estN=estN,pval=pval)
  }
  if(!printest) out <- pval
  return(out)
}


#' Power for Hypothesis Testing Using the Bailey Estimator
#' @description Approximates the power of a hypothesis test of the Bailey
#'   estimator by means of many simulated draws from a specified alternative distribution, conditioned on sample sizes.
#' @param nullN The abundance given by the null hypothesis
#' @param trueN The assumed abundance for the power calculation
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param alpha The alpha level for the test
#' @param nsim Number of simulated values to draw.  Defaults to 10000.
#' @param alternative Direction of the alternative hypothesis.  Allowed values
#'   are \code{"less"}, \code{"greater"}, or \code{"2-sided"}.  Defaults to \code{"less"}.
#' @return The approximate power of the specified hypothesis test, for the specified alternative value.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NBailey}, \link{vBailey}, \link{seBailey}, \link{rBailey},
#'   \link{pBailey}, \link{ciBailey}
#' @examples
#' powBailey(nullN=500, trueN=400, n1=100, n2=100, nsim=1000)
#'
#' Ntotry <- seq(from=250, to=450, by=25)
#' pows <- sapply(Ntotry, function(x)
#'   powBailey(nullN=500, trueN=x, n1=100, n2=100, nsim=1000))
#' plot(Ntotry, pows)
#' @export
powBailey <- function(nullN,trueN,n1,n2,alpha=.05,nsim=10000,alternative="less") {
  reject <- NA
  sim <- rBailey(nsim,trueN,n1,n2)
  for(i in 1:length(sim)) {
    reject[i] <- (pBailey(estN=sim[i],nullN=nullN,n1,n2,nsim=nsim,alternative=alternative)<alpha)
  }
  return(mean(reject))
}



#' Power for Hypothesis Testing Using the Chapman Estimator
#' @description Approximates the power of a hypothesis test of the Chapman
#'   estimator by means of many simulated draws from a specified alternative distribution, conditioned on sample sizes.
#' @param nullN The abundance given by the null hypothesis
#' @param trueN The assumed abundance for the power calculation
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param alpha The alpha level for the test
#' @param nsim Number of simulated values to draw.  Defaults to 10000.
#' @param alternative Direction of the alternative hypothesis.  Allowed values
#'   are \code{"less"}, \code{"greater"}, or \code{"2-sided"}.  Defaults to \code{"less"}.
#' @return The approximate power of the specified hypothesis test, for the specified alternative value.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{vChapman}, \link{seChapman}, \link{rChapman},
#'   \link{pChapman}, \link{ciChapman}
#' @examples
#' powChapman(nullN=500, trueN=400, n1=100, n2=100, nsim=1000)
#'
#' Ntotry <- seq(from=250, to=450, by=25)
#' pows <- sapply(Ntotry, function(x)
#'   powChapman(nullN=500, trueN=x, n1=100, n2=100, nsim=1000))
#' plot(Ntotry, pows)
#' @export
powChapman <- function(nullN,trueN,n1,n2,alpha=.05,nsim=10000,alternative="less") {
  reject <- NA
  sim <- rChapman(nsim,trueN,n1,n2)
  for(i in 1:length(sim)) {
    reject[i] <- (pChapman(estN=sim[i],nullN=nullN,n1,n2,nsim=nsim,alternative=alternative)<alpha)
  }
  return(mean(reject))
}



#' Power for Hypothesis Testing Using the Petersen Estimator
#' @description Approximates the power of a hypothesis test of the Petersen
#'   estimator by means of many simulated draws from a specified alternative distribution, conditioned on sample sizes.
#' @param nullN The abundance given by the null hypothesis
#' @param trueN The assumed abundance for the power calculation
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param alpha The alpha level for the test
#' @param nsim Number of simulated values to draw.  Defaults to 10000.
#' @param alternative Direction of the alternative hypothesis.  Allowed values
#'   are \code{"less"}, \code{"greater"}, or \code{"2-sided"}.  Defaults to \code{"less"}.
#' @return The approximate power of the specified hypothesis test, for the specified alternative value.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NPetersen}, \link{vPetersen}, \link{sePetersen}, \link{rPetersen},
#'   \link{pPetersen}, \link{ciPetersen}
#' @examples
#' powPetersen(nullN=500, trueN=400, n1=100, n2=100, nsim=1000)
#'
#' Ntotry <- seq(from=250, to=450, by=25)
#' pows <- sapply(Ntotry, function(x)
#'   powPetersen(nullN=500, trueN=x, n1=100, n2=100, nsim=1000))
#' plot(Ntotry, pows)
#' @export
powPetersen <- function(nullN,trueN,n1,n2,alpha=.05,nsim=10000,alternative="less") {
  reject <- NA
  sim <- rPetersen(nsim,trueN,n1,n2)
  for(i in 1:length(sim)) {
    reject[i] <- (pPetersen(estN=sim[i],nullN=nullN,n1,n2,nsim=nsim,alternative=alternative)<alpha)
  }
  return(mean(reject))
}


#' Plotting the Density of a Vector of Discrete Values
#' @description Plots the empirical density of a vector of discrete values, approximating the probability mass function (pmf).  This can be considered a more appropriate alternative to \code{plot(density(x))} in the case of a vector with a discrete (non-continuous) support, such as that calculated by an abundance estimator.
#' @param x The vector of values to plot
#' @param xlab The X-axis label for plotting
#' @param ylab The Y-axis label for plotting
#' @param ... Additional plotting arguments
#' @author Matt Tyers
#' @examples
#' draws <- rChapman(length=100000, N=500, n1=100, n2=100)
#' plotdiscdensity(draws)  #plots the density of a vector of discrete values
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @export
plotdiscdensity <- function(x,xlab="value",ylab="density",...) {
  uniquevals <- sort(unique(x))
  props <- NA
  for(i in 1:length(uniquevals)) {
    props[i] <- length(x[x==uniquevals[i]])/length(x)
  }
  plot(uniquevals,props,xlab=xlab,ylab=ylab,...=...)
  for(i in 1:length(uniquevals)) lines(rep(uniquevals[i],2),c(0,props[i]))
}
#discdensity(rPetersen(1000000,10000,1000,1000))

