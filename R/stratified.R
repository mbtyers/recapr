#' Stratified Abundance Estimator
#' @description Calculates the value of the stratified estimator for abundance
#'   in a mark-recapture experiment, from vectors of sample sizes and number of
#'   recaptures, with each element corresponding to each sampling stratum.
#' @param n1 Vector of individuals captured and marked in the first sample, from
#'   each stratum
#' @param n2 Vector of individuals captured and marked in the second sample,
#'   from each stratum
#' @param m2 Vector of marked individuals recaptured in the second sample, from
#'   each stratum
#' @param estimator The type of estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}.  Default to
#'   \code{"Chapman"}.
#' @return The value of the stratified estimator
#' @note It is possible that even the stratified estimate may be biased if
#'   capture probabilities differ greatly between strata.  However, the bias in
#'   the stratified estimator will be much less than an estimator calculated
#'   without stratification.
#' @author Matt Tyers
#' @seealso \link{strattest}, \link{rstrat}, \link{vstrat},  \link{sestrat}, \link{cistrat},
#'   \link{NChapman}, \link{NPetersen}, \link{NBailey}
#' @examples
#' Nstrat(n1=c(100,200), n2=c(100,500), m2=c(10,10))
#' @export
Nstrat <- function(n1,n2,m2,estimator="Chapman") {
  if(!(length(n1)==length(n2)) | !(length(n1)==length(m2))) stop("n1, n2, and m2 vectors must be of equal length")
  if(!estimator %in% c("Chapman","Bailey","Petersen")) stop("invalid estimator")
  if(estimator=="Chapman") est <- mapply(NChapman, n1=n1, n2=n2, m2=m2)
  if(estimator=="Bailey") est <- mapply(NBailey, n1=n1, n2=n2, m2=m2)
  if(estimator=="Petersen") est <- mapply(NPetersen, n1=n1, n2=n2, m2=m2)
  sum <- sum(est)
  return(sum)
}

#' Estimated Variance of Stratified Abundance Estimator
#' @description Calculates the estimated variance of the stratified estimator
#'   for abundance in a mark-recapture experiment, from vectors of sample sizes
#'   and number of recaptures, with each element corresponding to each sampling
#'   stratum.
#' @param n1 Vector of individuals captured and marked in the first sample, from
#'   each stratum
#' @param n2 Vector of individuals captured and marked in the second sample,
#'   from each stratum
#' @param m2 Vector of marked individuals recaptured in the second sample, from
#'   each stratum
#' @param estimator The type of estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}.  Default to
#'   \code{"Chapman"}.
#' @return The estimated variance of the stratified estimator
#' @note It is possible that even the stratified estimate of abundance may be
#'   biased if capture probabilities differ greatly between strata.  However,
#'   the bias in the stratified estimator will be much less than an estimator
#'   calculated without stratification.
#' @note This function makes the naive assumption of independence between
#'   strata.  Caution is therefore recommended.
#' @author Matt Tyers
#' @seealso \link{strattest}, \link{Nstrat}, \link{rstrat},  \link{sestrat}, \link{cistrat},
#'   \link{NChapman}, \link{NPetersen}, \link{NBailey}
#' @examples
#' vstrat(n1=c(100,200), n2=c(100,500), m2=c(10,10))
#' @export
vstrat <- function(n1,n2,m2,estimator="Chapman") {
  if(!(length(n1)==length(n2)) | !(length(n1)==length(m2))) stop("n1, n2, and m2 vectors must be of equal length")
  if(!estimator %in% c("Chapman","Bailey","Petersen")) stop("invalid estimator")
  if(estimator=="Chapman") est <- mapply(vChapman, n1=n1, n2=n2, m2=m2)
  if(estimator=="Bailey") est <- mapply(vBailey, n1=n1, n2=n2, m2=m2)
  if(estimator=="Petersen") est <- mapply(vPetersen, n1=n1, n2=n2, m2=m2)
  return(sum(est))
}

#vstrat(n1=c(100,200),n2=c(100,500),m2=c(10,10))
#vChapman(100,100,10)+vChapman(200,500,10)


#' Standard Error of Stratified Abundance Estimator
#' @description Calculates the standard error of the stratified estimator for
#'   abundance in a mark-recapture experiment, from vectors of sample sizes and
#'   number of recaptures, with each element corresponding to each sampling
#'   stratum.
#' @param n1 Vector of individuals captured and marked in the first sample, from
#'   each stratum
#' @param n2 Vector of individuals captured and marked in the second sample,
#'   from each stratum
#' @param m2 Vector of marked individuals recaptured in the second sample, from
#'   each stratum
#' @param estimator The type of estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}.  Default to
#'   \code{"Chapman"}.
#' @return The standard error of the stratified estimator
#' @note It is possible that even the stratified estimate of abundance may be
#'   biased if capture probabilities differ greatly between strata.  However,
#'   the bias in the stratified estimator will be much less than an estimator
#'   calculated without stratification.
#' @note This function makes the naive assumption of independence between
#'   strata.  Caution is therefore recommended.
#' @author Matt Tyers
#' @seealso \link{strattest}, \link{Nstrat}, \link{rstrat},  \link{vstrat}, \link{cistrat},
#'   \link{NChapman}, \link{NPetersen}, \link{NBailey}
#' @examples
#' sestrat(n1=c(100,200), n2=c(100,500), m2=c(10,10))
#' @export
sestrat <- function(n1,n2,m2,estimator="Chapman") {
  if(!(length(n1)==length(n2)) | !(length(n1)==length(m2))) stop("n1, n2, and m2 vectors must be of equal length")
  if(!estimator %in% c("Chapman","Bailey","Petersen")) stop("invalid estimator")
  sum <- sqrt(vstrat(n1=n1,n2=n2,m2=m2,estimator=estimator))
  return(sum)
}

#sestrat(n1=c(100,200),n2=c(100,500),m2=c(10,10))
#sqrt(vChapman(100,100,10)+vChapman(200,500,10))


#' Random Draws from the Stratified Estimator
#' @description Returns a vector of random draws from the stratified estimator in a
#'   mark-recapture experiment, given values of the true abundance and the
#'   sample size in both events.  The function first simulates a vector of
#'   recaptures (m2) for each stratum, and then uses these to
#'   compute a vector of draws from the estimator.
#'
#'   It may prove useful to investigate the behavior of the stratified estimator under relevant scenarios.
#'
#' If capture probabilities (\code{p1} and/or \code{p2}) are specified instead of sample size(s), the sample size(s) will first be drawn from a binomial distribution, then the number of recaptures.  If both sample size and capture probability are specified for a given sampling event, only the sample size will be used.
#' @param length The length of the random vector to return.
#' @param N A vector of values of the true abundance for each stratum.
#' @param n1 A vector of the number of individuals captured and marked in the first sample, for each stratum.
#' @param n2 A vector of the number of individuals captured in the second sample, for each stratum.
#' @param p1 Alternately, a vector of probabilities of capture for the first event, for each stratum.
#' @param p2 Alternately, a vector of probabilities of capture for the second event, for each stratum.
#' @param estimator The type of estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}.  Default to
#'   \code{"Chapman"}.
#' @return A vector of random draws from the stratified estimator
#' @author Matt Tyers
#' @seealso \link{strattest}, \link{Nstrat},  \link{vstrat}, \link{cistrat},\link{NChapman}, \link{NPetersen}, \link{NBailey}
#' @importFrom stats rhyper
#' @examples
#' draws <- rstrat(length=100000, N=c(5000,10000), n1=c(500,200), n2=c(500,200))
#' plotdiscdensity(draws)  #plots the density of a vector of discrete values
#' mean(draws)
#' @export
rstrat <- function(length,N,n1=NULL,n2=NULL,p1=NULL,p2=NULL,estimator="Chapman") {
  if(!is.null(n1) & !is.null(p1)) warning("both n1 and p1 specified - only n1 used")
  if(!is.null(n2) & !is.null(p2)) warning("both n2 and p2 specified - only n2 used")
  if(!is.null(n1) & !is.null(n2)) {
    if(!(length(n1)==length(n2)) | !(length(n1)==length(N))) stop("n1, n2, and N vectors must be of equal length")
  }
  if(!estimator %in% c("Chapman","Bailey","Petersen")) stop("invalid estimator")
  Nhat <- 0
  for(i in 1:length(N)) {
    if(!is.null(n1)) n1fill <- n1[i]
    if(!is.null(n2)) n2fill <- n2[i]
    if(is.null(n1) & is.null(n2)) {
      if(is.null(p1)) stop("either n1 or p1 vectors must be specified")
      if(is.null(p2)) stop("either n2 or p2 vectors must be specified")
    }
    if(is.null(n1)) {
      if(is.null(p1)) stop("either n1 or p1 vectors must be specified")
      if(!is.null(n2)) {
        if(!(length(p1)==length(n2)) | !(length(p1)==length(N))) stop("p1, n2, and N vectors must be of equal length")
      }
      if(is.null(n2)) {
        if(!(length(p1)==length(p2)) | !(length(p1)==length(N))) stop("p1, p2, and N vectors must be of equal length")
      }
      n1fill <- rbinom(length,size=N[i],prob=p1[i])
    }
    if(is.null(n2)) {
      if(is.null(p2)) stop("either n2 or p2 vectors must be specified")
      if(!is.null(n1)) {
        if(!(length(p2)==length(n1)) | !(length(p2)==length(N))) stop("p2, n1, and N vectors must be of equal length")
      }
      if(is.null(n1)) {
        if(!(length(p2)==length(p1)) | !(length(p2)==length(N))) stop("p1, p2, and N vectors must be of equal length")
      }
      n2fill <- rbinom(length,size=N[i],prob=p2[i])
    }
    if(estimator=="Chapman") {
      m2 <- rhyper(length, n1fill, N[i]-n1fill, n2fill)
      Nhat <- Nhat + NChapman(n1=n1fill, n2=n2fill, m2=m2)
    }
    if(estimator=="Petersen") {
      m2 <- rhyper(length, n1fill, N[i]-n1fill, n2fill)
      Nhat <- Nhat + NPetersen(n1=n1fill, n2=n2fill, m2=m2)
    }
    if(estimator=="Bailey") {
      m2 <- rbinom(length, n2fill, n1fill/N[i])
      Nhat <- Nhat + NBailey(n1=n1fill, n2=n2fill, m2=m2)
    }
  }
  return(Nhat)
  # if(!(length(n1)==length(n2)) | !(length(n1)==length(N))) stop("n1, n2, and N vectors must be of equal length")
  # if(!estimator %in% c("Chapman","Bailey","Petersen")) stop("invalid estimator")
  # Nhat <- 0
  # for(i in 1:length(N)) {
  #   if(estimator=="Chapman") {
  #     m2 <- rhyper(length, n1[i], N[i]-n1[i], n2[i])
  #     Nhat <- Nhat + NChapman(n1=n1[i], n2=n2[i], m2=m2)
  #   }
  #   if(estimator=="Petersen") {
  #     m2 <- rhyper(length, n1[i], N[i]-n1[i], n2[i])
  #     Nhat <- Nhat + NPetersen(n1=n1[i], n2=n2[i], m2=m2)
  #   }
  #   if(estimator=="Bailey") {
  #     m2 <- rbinom(length, n2[i], n1[i]/N[i])
  #     Nhat <- Nhat + NBailey(n1=n1[i], n2=n2[i], m2=m2)
  #   }
  # }
  # return(Nhat)
}




#' Confidence Intervals for the Stratified Estimator
#' @description Calculates approximate confidence intervals(s) for the
#'   Stratified estimator, using bootstrapping, the Normal approximation, or
#'   both.
#'
#'   The bootstrap interval is created by resampling the data in the second
#'   sampling event, with replacement for each stratum; that is, drawing
#'   bootstrap values of m2 from a binomial distribution with probability
#'   parameter m2/n2.
#' @param n1 Number of individuals captured and marked in the first sample
#' @param n2 Number of individuals captured in the second sample
#' @param m2 Number of marked individuals recaptured in the second sample
#' @param conf The confidence level of the desired intervals.  Defaults to 0.95.
#' @param method Which method of confidence interval to return.  Allowed values
#'   are \code{"norm"}, \code{"boot"}, or \code{"both"}.  Defaults to
#'   \code{"both"}.
#' @param bootreps Number of bootstrap replicates to use.  Defaults to 10000.
#' @param estimator The type of estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}.  Default to
#'   \code{"Chapman"}.
#' @param useChapvar Whether to use the Chapman estimator variance instead of
#'   the Petersen estimator variance for the normal-distribution interval, if
#'   \code{"method"} is set to  \code{"Petersen"}. Defaults to \code{FALSE}.
#' @return A list with the abundance estimate and confidence interval bounds for
#'   the normal-distribution and/or bootstrap confidence intervals.
#' @note Both the bootstrap and the normal approximation intervals make the
#'   naive assumption of independence between strata, which may not be the case.
#'   The user therefore cautioned, and is encouraged to investigate the coverage
#'   under relevant scenarios.
#' @author Matt Tyers
#' @seealso \\link{strattest}, \link{Nstrat}, \link{rstrat},  \link{vstrat}, \link{sestrat},
#'   \link{NChapman}, \link{NPetersen}, \link{NBailey}
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats quantile
#' @examples
#' cistrat(n1=c(100,200), n2=c(100,500), m2=c(10,10))
#' @export
cistrat <- function(n1,n2,m2,conf=0.95, method="both", bootreps=10000, estimator="Chapman", useChapvar=FALSE) {
  if(!(length(n1)==length(n2)) | !(length(n1)==length(m2))) stop("n1, n2, and m2 vectors must be of equal length")
  if(!estimator %in% c("Chapman","Bailey","Petersen")) stop("invalid estimator")
  if(!any(method==c("boot","norm","both"))) stop("invalid method specified")
  bounds <- c((1-conf)/2,1-((1-conf)/2))
  out <- list()

  out$Nstrat <- Nstrat(n1=n1, n2=n2, m2=m2, estimator=estimator)
  Nhat_by_strat <- NA
  for(i in 1:length(n1)) {
    if(estimator=="Chapman") Nhat_by_strat[i] <- NChapman(n1=n1[i],n2=n2[i],m2=m2[i])
    if(estimator=="Petersen") Nhat_by_strat[i] <- NPetersen(n1=n1[i],n2=n2[i],m2=m2[i])
    if(estimator=="Bailey") Nhat_by_strat[i] <- NBailey(n1=n1[i],n2=n2[i],m2=m2[i])
  }
  out$Nhat_by_strat <- Nhat_by_strat

  if(method=="norm" | method=="both") {
    se <- ifelse(useChapvar,sestrat(n1=n1,n2=n2,m2=m2,estimator="Chapman"),sestrat(n1=n1,n2=n2,m2=m2,estimator=estimator))
    ciNorm <- out$Nstrat + qnorm(bounds)*se
    out$ciNorm_strat <- ciNorm
  }

  if(method=="boot" | method=="both") {
    Nhatboot <- 0
    for(i in 1:length(n1)) {
      m2boot <- rbinom(bootreps,size=n2[i],prob=(m2[i]/n2[i]))
      if(estimator=="Bailey") Nhatboot <- Nhatboot + NBailey(n1=n1[i],n2=n2[i],m2=m2boot)
      if(estimator=="Chapman") Nhatboot <- Nhatboot + NChapman(n1=n1[i],n2=n2[i],m2=m2boot)
      if(estimator=="Petersen") Nhatboot <- Nhatboot + NPetersen(n1=n1[i],n2=n2[i],m2=m2boot)
    }
    #Nhatboot <- Nhatboot +(out$Nstrat-mean(Nhatboot))    # trying to bootstrap-unbias it
    out$ciBoot_strat <- unname(quantile(Nhatboot,bounds,na.rm=TRUE))
    #out$bootmean <- mean(Nhatboot)
  }
  return(out)
}

#cistrat(n1=c(100,200),n2=c(100,500),m2=c(10,10))
#ciChapman(100,100,10)
#ciChapman(200,500,10)
#ciChapman(200,700,20)

#
# N <- c(10000,2000,10000)
# n1 <- c(1000,100,100)
# n2 <- c(1000,100,100)
# nsim <- 10000
# m21 <- rhyper(nsim,n1[1],N[1]-n1[1],n2[1])
# m22 <- rhyper(nsim,n1[2],N[2]-n1[2],n2[2])
# m23 <- rhyper(nsim,n1[3],N[3]-n1[3],n2[3])
#
# # less <- more <- gotit <- rep(F,nsim)
# # bootmean <- NA
# # for(i in 1:nsim) {
# #   x <- cistrat(n1=n1,n2=n2,m2=c(m21[i],m22[i],m23[i]),conf=.9)
# #   ci <- x$ciBoot_strat
# #   if(sum(N)<ci[1]) less[i] <- T
# #   if(sum(N)>ci[2]) more[i] <- T
# #   if(sum(N)>=ci[1] & sum(N)<=ci[2]) gotit[i] <- T
# #   bootmean[i] <- x$bootmean
# # }
# # mean(less)
# # mean(gotit)
# # mean(more)
# # mean(bootmean)
# # #
# # # Nstrat <- NChapman(n1[1],n2[1],m21) + NChapman(n1[2],n2[2],m22) + NChapman(n1[3],n2[3],m23)
# # # mean(Nstrat)
# # # sum(N)
# #
#
# sum(N)
# Nhat <- Nhat_unb <- Nhathat<- NA
# for(i in 1:nsim) {
#   asdf <- Nstrat(n1=n1,n2=n2,m2=c(m21[i],m22[i],m23[i]))
#   Nhat[i] <- asdf$Nhat
#   Nhat_unb[i] <- asdf$Nhat_unb
#   Nhathat[i] <- NChapman(n1=sum(n1),n2=sum(n2),m2=sum(c(m21[i],m22[i],m23[i])))
#   if(20*i/nsim==floor(20*i/nsim)) cat(100*i/nsim,"% ... ")
# }
# mean(Nhathat)
# mean(Nhat)
# mean(Nhat_unb)
# hist(Nhathat)
# hist(Nhat)
# hist(Nhat_unb)
# sd(Nhathat)
# sd(Nhat)
# sd(Nhat_unb)
#
#
# nsim <- 1000
# N <- c(10000,15000)
# p1 <- p2 <- seq(.1,.9,by=.1)
# n11 <- n12 <- p1*N[1] #event, population
# n21 <- n22 <- .5*p2*N[2]
# Bhat <- Bhat_mean <- Bhat_median <- Bhathat <- matrix(NA,nrow=length(p1),ncol=length(p2))
# for(i in 1:length(p1)) {
#   for(j in 1:length(p2)) {
#     m21 <- rhyper(nsim,n11[i],N[1]-n11[i],n21[j])
#     m22 <- rhyper(nsim,n12[i],N[2]-n12[i],n22[j])
#     Nhat <- Nhat_mean <- Nhat_median <- Nhathat <- NA
#     for(k in 1:nsim) {
#       asdf <- Nstrat(n1=c(n11[i],n12[i]),n2=c(n21[j],n22[j]),m2=c(m21[k],m22[k]))
#       Nhat[k] <- asdf$Nhat
#       Nhat_mean[k] <- asdf$Nhat_unb
#       Nhat_median[k] <- asdf$Nhat_medunb
#       Nhathat[k] <- NChapman(n1=sum(n11[i],n12[i]),n2=sum(n21[j],n22[j]),m2=sum(c(m21[k],m22[k])))
#     }
#     Bhat[i,j] <- mean(Nhat,na.rm=TRUE)-sum(N)
#     Bhat_mean[i,j] <- mean(Nhat_mean,na.rm=TRUE)-sum(N)
#     Bhat_median[i,j] <- mean(Nhat_median,na.rm=TRUE)-sum(N)
#     Bhathat[i,j] <- mean(Nhathat,na.rm=TRUE)-sum(N)
#   }
#   print(i)
# }
# Bhat
# Bhat_mean
# Bhat_median
# Bhathat
#
# image(Bhat)
# plot(as.vector(Bhat),as.vector(Bhathat))
#
# plot(as.vector(Bhat),as.vector(Bhat_mean))
# abline(c(1,1))
#
# nsim <- 1000
# N1 <- seq(1000,10000,by=1000)#rep(10000,10)
# N2 <- seq(1000,10000,by=1000)#rep(10000,10) #seq(1000,10000,by=1000) #c(10000,15000)
# p1 <- seq(.05,.5,by=.05)
# p2 <- rep(.1,10)#seq(.05,.5,by=.05)#.1
# n11 <- n12 <- p1*N1 #population,event
# n21 <- n22 <- p2*N2
# Bhat <- Bhat_mean <- Bhat_median <- Bhathat <- matrix(NA,nrow=length(N1),ncol=length(N2))
# for(i in 1:length(N1)) {
#   for(j in 1:length(N2)) {
#     m12 <- rhyper(nsim,n11[i],N1[i]-n11[i],n12[i])
#     m22 <- rhyper(nsim,n21[j],N2[j]-n21[j],n22[j])
#     Nhat <- Nhat_mean <- Nhat_median <- Nhathat <- NA
#     for(k in 1:nsim) {
#       asdf <- Nstrat(n1=c(n11[i],n21[j]),n2=c(n12[i],n22[j]),m2=c(m12[k],m22[k]))
#       Nhat[k] <- asdf$Nhat
#       Nhat_mean[k] <- asdf$Nhat_unb
#       Nhat_median[k] <- asdf$Nhat_medunb
#       Nhathat[k] <- NChapman(n1=sum(n11[i],n21[j]),n2=sum(n12[i],n22[j]),m2=sum(c(m12[k],m22[k])))
#     }
#     Bhat[i,j] <- mean(Nhat,na.rm=TRUE)-(N1[i]+N2[j])#sum(N)
#     Bhat_mean[i,j] <- mean(Nhat_mean,na.rm=TRUE)-(N1[i]+N2[j])#sum(N)
#     Bhat_median[i,j] <- mean(Nhat_median,na.rm=TRUE)-(N1[i]+N2[j])#sum(N)
#     Bhathat[i,j] <- mean(Nhathat,na.rm=TRUE)-(N1[i]+N2[j])#sum(N)
#   }
#   print(i)
# }
# Bhat
# Bhat_mean
# Bhat_median
# Bhathat
#
# image(Bhat)
# plot(as.vector(Bhat),as.vector(Bhathat))
#
# plot(as.vector(Bhat),as.vector(Bhat_median))
# abline(c(1,1))
# plot(as.vector(Bhat),as.vector(Bhat_mean))
# abline(c(1,1))
#
#
#
# N1 <- 10000
# N2 <- 10000
# n11 <- n12 <- 1000
# n21 <- n22 <- 500
# m12 <- rhyper(nsim,n11,N1-n11,n12)
# m22 <- rhyper(nsim,n21,N2-n21,n22)
# Nhat <- Nhat_mean <- Nhat_median <- Nhathat <- NA
# for(k in 1:nsim) {
#   asdf <- Nstrat(n1=c(n11,n21),n2=c(n12,n22),m2=c(m12[k],m22[k]))
#   Nhat[k] <- asdf$Nhat
#   Nhat_mean[k] <- asdf$Nhat_unb
#   Nhat_median[k] <- asdf$Nhat_medunb
#   Nhathat[k] <- NChapman(n1=sum(n11,n21),n2=sum(n12,n22),m2=sum(c(m21[k],m22[k])))
#   if(20*k/nsim==floor(20*k/nsim)) cat(100*k/nsim,"% ... ")
# }
# mean(Nhat,na.rm=TRUE)-(N1+N2)#(N1[i]+N2[j])#sum(N)
# mean(Nhat_mean,na.rm=TRUE)-(N1+N2)#(N1[i]+N2[j])#sum(N)
# mean(Nhat_median,na.rm=TRUE)-(N1+N2)#(N1[i]+N2[j])#sum(N)
# mean(Nhathat,na.rm=TRUE)-(N1+N2)#(N1[i]+N2[j])#sum(N)
