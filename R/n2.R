#' Mark-Recapture Sample Size, Robson-Regier
#' @description Calculates minimum sample size for one sampling event in a
#'   Petersen mark-recapture experiment, given the sample size in the other
#'   event and an best guess at true abundance.
#' @param N The best guess at true abundance
#' @param n1 The size of the first (or second) sampling event
#' @param conf A vector of the desired levels of confidence to investigate.
#'   Allowed values are any of \code{c(0.99,0.95,0.9,0.85,0.8,0.75)}.  Defaults to
#'   all of \code{c(0.99,0.95,0.85,0.8,0.75)}.
#' @param acc A vector of the desired levels of relative accuracy to
#'   investigate.  Allowed values are any of
#'   \code{c(0.5,0.25,0.2,0.15,0.1,0.05,0.01)}.  Defaults to all of
#'   \code{c(0.5,0.25,0.2,0.15,0.1,0.05,0.01)}.
#' @return A list of minimum sample sizes.  Each list element corresponds to a
#'   unique level of confidence, and is defined as a data frame with each row
#'   corresponding to a unique value of relative accuracy.  Two minimum sample
#'   sizes are given: one calculated from the sample size provided for the other
#'   event, and the other calculated under n1=n2, the most efficient scenario.
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @note It is possible that the sample size - accuracy relationship will be better illustrated using \link{plotn2sim}.
#' @references Robson, D. S., and H. A. Regier.  1964.  Sample size in Petersen
#'   mark-recapture experiments.  Transactions of the American FisheriesSociety
#'   93:215-226.
#' @author Matt Tyers
#' @seealso \link{plotn2sim}, \link{plotn1n2simmatrix}
#' @importFrom stats rhyper
#' @examples
#' n2RR(N=1000, n1=100)
#' @export
n2RR <- function(N, n1, conf=c(0.99,0.95,0.9,0.85,0.8,0.75),acc=c(0.5,0.25,0.2,0.15,0.1,0.05,0.01)) {
  out <- list()
  accall <- c(0.5,0.25,0.2,0.15,0.1,0.05,0.01)
  conf <- conf[conf %in% c(0.99,0.95,0.9,0.85,0.8,0.75)]
  acc <- acc[acc %in% accall]
  if(length(acc)==0) stop("invalid acc")
  if(length(conf)==0) stop("invalid conf")
  acc <- sort(acc,decreasing=T)
  for(i in 1:length(conf)) {
    Draw <- (conf[i]==0.99)*c(48.707,135.529229,196.2951017,325.8179089,694.4679116,2684.858928,66380.75202) +
      (conf[i]==0.95)*c(24.35001031,69.83385388,103.9898307,178.3155987,391.4504252,1543.721381,38421.81) +
      (conf[i]==0.9)*c(14.7893,45.79545303,69.92066847,122.361479,272.5531006,1084.149431,27057.47792) +
      (conf[i]==0.85)*c(9.74015,33.5685,52.11740329,92.33768116,207.4134585,829.0626252,20722.71022) +
      (conf[i]==0.8)*c(6.64835,25.80786976,40.54417323,72.44665423,163.6673542,656.3687191,16423.22022) +
      (conf[i]==0.75)*c(4.700876432,20.3312336,32.2203623,57.93630836,131.443365,528.4288547,13232.17423)
    D <- Draw[accall %in% acc]
    n2_1 <- ceiling(N/(n1*(N-1)/((N-n1)*D)+1))
    n2_2 <- ceiling(N/(sqrt((N-1)/D)+1))
    df <- data.frame(acc,n2_1,n2_2)
    dimnames(df)[[2]] <- c("Rel acc","n2 from specified n1","n if n1=n2")
    out[[i]] <- df
    names(out)[[i]] <- paste0("conf_",conf[i])
  }
  return(out)
}

#n2RR(1000,100,conf="frogs",acc=c(.2,.5,.333,"poo"))

#' Mark-Recapture Sample Size Via Simulation
#' @description Produces a plot of the simulated relative accuracy of a
#'   mark-recapture abundance estimator for various sample sizes.  This may be a
#'   better representation of the sample size - accuracy relationship than that
#'   provided by \link{n2RR}.
#' @param N The best guess at true abundance
#' @param n1 The size of the first (or second) sampling event
#' @param conf A vector of the desired levels of confidence to investigate.
#'   Allowed values are any of \code{c(0.99,0.95,0.85,0.8,0.75)}.  Defaults to
#'   all of \code{c(0.99,0.95,0.85,0.8,0.75)}.
#' @param n2range A two-element vector describing the range of sample sizes to
#'   investigate.  If the default (\code{NULL}) is accepted, an appropriate
#'   value will be chosen.
#' @param n2step The step size between sample sizes to investigate.  If the
#'   default (\code{NULL}) is accepted, an appropriate value will be chosen.
#' @param estimator The abundance estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}. Defaults to
#'   \code{"Chapman"}.
#' @param nsim The number of replicates.  Defaults to 10000.
#' @param accrange The maximum level of relative accuracy for plotting.
#'   Defaults to 1.
#' @param ... Additional plotting parameters
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{n2RR}, \link{plotn1n2simmatrix}
#' @importFrom stats quantile
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @importFrom graphics grid
#' @importFrom graphics par
#' @examples
#' plotn2sim(N=1000, n1=100)
#' @export
plotn2sim <- function(N, n1, conf=c(0.99,0.95,0.85,0.8,0.75),n2range=NULL,n2step=NULL,estimator="Chapman",nsim=10000,accrange=1,...) {
  if(is.null(n2range)) n2range <- c(0,2*as.numeric(n2RR(N=N,n1=n1,conf=.95,acc=.1)[[1]])[3])
  if(is.null(n2step)) n2step <- max(1,floor((n2range[2]-n2range[1])/200))
  if(length(n2range)!=2) stop("n2range must have two elements (min and max n2 to plot)")
  whichn2 <- seq(n2range[1],n2range[2],by=n2step)
  if(sum(conf %in% c(0.99,0.95,0.85,0.8,0.75)) < 1) stop("invalid conf")
  conftry <- sort(conf[conf %in% c(0.99,0.95,0.85,0.8,0.75)],decreasing=T)
  accsimcalc <- matrix(NA,nrow=length(whichn2),ncol=(2*length(conftry)+1))
  accsimcalc[,1] <- whichn2
  i <- 1
  for(nn22 in whichn2) {
    if(estimator=="Chapman") N.Chap <- rChapman(length=nsim, N=N, n1=n1, n2=nn22)
    if(estimator=="Petersen") N.Chap <- rPetersen(length=nsim, N=N, n1=n1, n2=nn22)
    if(estimator=="Bailey") N.Chap <- rBailey(length=nsim, N=N, n1=n1, n2=nn22)
    # m2 <- rhyper(input$nsim,input$n1_a,input$N_a-input$n1_a,nn22)
    # N.Chap <- (input$n1_a+1)*(nn22+1)/(m2+1) - 1
    for(j in 2:(length(conftry)+1)) {
      quantilesN <- quantile(N.Chap,c((1-conftry[j-1])/2,(1+conftry[j-1])/2))
      quantilesRP <- (quantilesN-N)/N
      accsimcalc[i,j] <- quantilesRP[1]
      accsimcalc[i,(j+length(conftry))] <- quantilesRP[2]
    }
    i <- i+1
  }
  # accsimcalc
  cols <- c("#CC0000FF", "#CCCC00FF", "#00CC00FF", "#00CCCCFF", "#0000CCFF", "#CC00CCFF")
  plot(NA,xlim=c(accsimcalc[1,1],accsimcalc[1,1]+1*(accsimcalc[dim(accsimcalc)[1],1]-accsimcalc[1,1])),
       ylim=c(-(accrange),(accrange)),xlab="sample size n2",ylab="Simulated relative accuracy",
       main=c(paste0("Anticipated N: ",N,paste(", Selected n1:",n1))),...=...)
  grid()
  for(i in 2:(length(conftry)+1)) {
    lines(accsimcalc[,1],accsimcalc[,i],col=cols[i-1],lwd=2)
    lines(accsimcalc[,1],accsimcalc[,(i+length(conftry))],col=cols[i-1],lwd=2)
  }
  legend(par("usr")[1]+0.82*(par("usr")[2]-par("usr")[1]),
         par("usr")[4],
         legend=conftry,lwd=2,col=cols,title="conf (1-alpha)")


}

#n2sim(N=1000,n1=100)

#' Mark-Recapture Sample Size Via Sim, Both Samples
#' @description Produces a plot of the simulated relative accuracy of a
#'   mark-recapture abundance estimator for various sample sizes in both events.  This may be a
#'   useful exploration, but it is possible that \link{plotn2sim} may be more informative.
#' @param N The best guess at true abundance
#' @param conf The desired level of confidence to investigate.
#'   Defaults to 0.95.
#' @param nrange A two-element vector describing the range of sample sizes to
#'   investigate.  If the default (\code{NULL}) is accepted, an appropriate
#'   value will be chosen.
#' @param nstep The step size between sample sizes to investigate.  If the
#'   default (\code{NULL}) is accepted, an appropriate value will be chosen.
#' @param estimator The abundance estimator to use.  Allowed values are
#'   \code{"Chapman"}, \code{"Petersen"}, and \code{"Bailey"}. Defaults to
#'   \code{"Chapman"}.
#' @param nsim The number of replicates.  Defaults to 10000.
#' @param ... Additional plotting parameters
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{n2RR}, \link{plotn2sim}
#' @importFrom stats quantile
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics legend
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices rainbow
#' @importFrom graphics grid
#' @examples
#' plotn1n2simmatrix(N=10000, nsim=2000)
#' @export
plotn1n2simmatrix <- function(N,conf=0.95,nrange=NULL,nstep=NULL,estimator="Chapman",nsim=10000,...) {
  if(is.null(nrange)) nrange <- c(0,2*as.numeric(n2RR(N=N,n1=100,conf=conf,acc=.1)[[1]])[3])
  if(length(nrange)!=2) stop("nrange must have two elements (min and max n2 to plot)")
  if(is.null(nstep)) nstep <- max(1,floor((nrange[2]-nrange[1])/50))
  n1 <- n2 <- seq(from=nrange[1],to=nrange[2],by=nstep)
  simacc <- matrix(NA,nrow=length(n1),ncol=length(n2))
  for(i in 1:length(n1)) {
    for(j in 1:length(n2)) {
      if(estimator=="Chapman") Nhat <- rChapman(length=nsim,N=N,n1=n1[i],n2=n2[j])
      if(estimator=="Petersen") Nhat <- rPetersen(length=nsim,N=N,n1=n1[i],n2=n2[j])
      if(estimator=="Bailey") Nhat <- rBailey(length=nsim,N=N,n1=n1[i],n2=n2[j])
      quantilesN <- quantile(Nhat,c((1-conf)/2,(1+conf)/2),na.rm=T)
      quantilesRP <- (quantilesN-N)/N
      simacc[i,j] <- max(abs(quantilesRP))
    }
  }
  #contour(simacc,levels=c(0.5,0.25,0.2,0.15,0.1,0.05,0.01))
  #cols <- c("#CC0000FF", "#CCCC00FF", "#00CC00FF", "#00CCCCFF", "#0000CCFF", "#CC00CCFF")
  accbreaks <- c(10,0.5,0.25,0.2,0.15,0.1,0.05,0.01,0)
  cols <- adjustcolor(rainbow(length(accbreaks)+1),alpha.f=.85,red.f=.85,green.f=.85,blue.f=.85)
  colmat <- cols[as.numeric(cut(simacc,accbreaks))]
  y <- rep(n2,length(n1))
  x <- NULL
  for(i in 1:length(n2)) x <- c(x,rep(n1[i],length(n2)))
  xmat <- matrix(x,nrow=length(n1),ncol=length(n2))
  ymat <- matrix(y,nrow=length(n1),ncol=length(n2))
  #plot(xmat,ymat,col=colmat)
  plot(n1,n2,col="white",...=...)
  xleft<-xmat-nstep/2
  xright<-xmat+nstep/2
  ybottom <- ymat - nstep/2
  ytop <- ymat + nstep/2
  rect(xleft,ybottom,xright,ytop,col=colmat,border=NA)
  grid(col=1)
  #throwaway <- c(.003,.03,.07,.12,.17,.22,.3,.7)
  legend(0,max(n2),legend=c("0.00 - 0.01", "0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.15", "0.15 - 0.20", "0.20 - 0.25", "0.25 - 0.50", "0.50 +"), col=cols, pch=15,title=paste0("rel acc, conf=",conf))
}

#n12simmatrix(N=10000,nsim=20000,nrange=c(900,1100))
