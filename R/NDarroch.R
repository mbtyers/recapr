#' Spatially or Temporally Stratified Abundance Est (Darroch)
#' @description Computes abundance estimates and associated variance in the event of spatial or temporal stratification, or in any stratification in which individuals can move between strata.  Marking (event 1) and recapture (event 2) strata do not need to be the same.
#'
#' Inputs are vectors of total event 1 and 2 sample sizes, and either vectors of event 1 and 2 strata corresponding to each recaptured individual, or a matrix of total number of recaptures for each combination of event 1 and event 2 strata.
#'
#' Implementation is currently using Darroch's method, and will only accept non-singular input matrices.
#' @references Darroch, J.N. (1961). The two-sample capture-recapture census when tagging and sampling are stratified.  \emph{Biometrika} \strong{48}, 241-60.
#' @return A numeric list, giving the strata matrix if originally given in vector form, abundance estimates and standard errors by event 1 and event 2 strata, and the total abundance estimate and standard error.
#' @param n1counts A vector of the total sample sizes in the first event, by
#'   strata.  For example, setting \code{n1counts=c(20,30,40)} would imply 20
#'   individuals captured and marked in stratum 1, 30 in stratum 2, and 40 in
#'   stratum 3.
#' @param n2counts A vector of the total sample sizes in the second event, by
#'   strata.
#' @param m2strata1 A vector of the first-event stratum membership of each
#'   recaptured individual.  Only values \code{1, 2, 3, ...} are allowed.  May
#'   be used together with \code{m2strata2} instead of \code{stratamat}.
#' @param m2strata2 A vector of the second-event stratum membership of each
#'   recaptured individual.  Only values \code{1, 2, 3, ...} are allowed.  May
#'   be used together with \code{m2strata1} instead of \code{stratamat}.
#' @param stratamat A matrix specifying the number of recaptures in each
#'   combination of event 1 and event 2 strata, with rows corresponding to event
#'   1 strata and columns corresponding to event 2 strata.  May be used instead
#'   of \code{m2strata1} and \code{m2strata2}.
#' @author Matt Tyers
#' @seealso \link{consistencytest}
#' @importFrom MASS ginv
#' @examples
#' mat <- matrix(c(59,30,1,45,280,38,0,42,25), nrow=3, ncol=3, byrow=TRUE)
#'
#' NDarroch(n1counts=c(484,1468,399), n2counts=c(847,6616,2489), stratamat=mat)
#' @export
NDarroch <- function(n1counts, n2counts, m2strata1=NULL, m2strata2=NULL, stratamat=NULL) {
  out <- list()
  if(!is.null(stratamat)) {
    if((nrow(stratamat) != length(n1counts)) | (ncol(stratamat) != length(n2counts))) stop("Dimension mismatch in stratamat - must have rows corresponding to n1 strata and columns corresponding to n1 strata")
  }
  if(is.null(stratamat)) {
    if(is.null(m2strata1) | is.null(m2strata2)) stop("recapture strata must be specified, either with stratamat or m2strata1 and m2strata2 together")
    if(any(!is.numeric(m2strata1))) stop("m2strata1 values must be positive integers")
    if(any(!is.numeric(m2strata2))) stop("m2strata2 values must be positive integers")
    if(sum(m2strata1-round(m2strata1))!=0) stop("m2strata1 values must be positive integers")
    if(sum(m2strata2-round(m2strata2))!=0) stop("m2strata2 values must be positive integers")
    if(any(m2strata1<1)) stop("m2strata1 values must be positive integers")
    if(any(m2strata2<1)) stop("m2strata2 values must be positive integers")
    if(max(m2strata1)>length(n1counts)) stop("strata values larger than count vector in event 1")
    if(max(m2strata2)>length(n2counts)) stop("strata values larger than count vector in event 2")
    if(length(m2strata1)!=length(m2strata2)) stop("m2strata1 and m2strata2 must be the same length (each element corresponding to an individual)")
    stratamat <- matrix(NA,nrow=length(n1counts),ncol=length(n2counts))
    for(i in 1:length(n1counts)) {
      for(j in 1:length(n2counts)) {
        stratamat[i,j] <- sum(m2strata1==i & m2strata2==j)
      }
    }
    out$stratamat <- stratamat
  }
  mjdot <- colSums(stratamat)
  midot <- rowSums(stratamat)
  f <- function(m) class(try(solve(m),silent=T))=="matrix"
  if(!f(stratamat)) stop("movement matrix is non-invertable.") ###
  Minv <- ginv(stratamat)
  a <- n1counts
  u <- n2counts
  M <- stratamat
  rr <- Minv %*% a
  p <- 1/rr
  D_a <- diag(a)
  D_r <- diag(as.vector(rr))
  D_u <- diag(u)
  if(!f(D_a)) stop("D_a matrix is non-invertable.") ###
  D_ainv <- ginv(D_a)
  thetahat <- D_ainv %*% M %*% D_r
  if(!f(D_a)) stop("theta_hat matrix is non-invertable.") ###
  thetahatinv <- ginv(thetahat)
  Uhat <- D_u %*% rr
  Nhat <- Uhat + mjdot*rr
  Nhat_sum <- sum(Nhat)

  Uhatstar <- t(thetahatinv) %*% Uhat
  Nhatstar <- Uhatstar + a
  Nhatstar_sum <- sum(Nhatstar)

  mu <- thetahat %*% rr -1
  D_mu <- diag(as.vector(mu))

  X <- thetahatinv %*% D_mu %*% D_ainv %*% t(thetahatinv)
  Sig_r <- D_r %*% X %*% D_r
  Ident <- diag(nrow=nrow(D_r))
  SigUhat <- (D_u %*% Sig_r %*% D_u) + (D_u %*% D_r %*% (D_r-Ident))
  ones <- rep(1,nrow(D_r))
  varNhat_sum <- t(ones) %*% SigUhat %*% ones

  varNhat <- diag(SigUhat)

  B <- matrix(NA,nrow=ncol(thetahat),ncol=ncol(thetahat))
  D_A <- matrix(0,nrow=ncol(thetahat),ncol=ncol(thetahat))
  for(j in 1:ncol(thetahat)) {
    for(k in 1:ncol(thetahat)) {
      B[j,k] <- sum(thetahat[,j]*thetahat[,k]*(Uhatstar^2)/a)
    }
    D_A[j,j] <- sum(thetahat[,j]*(as.vector(Uhatstar)^2)/a)
  }
  D_p <- diag(as.vector(p))
  Ident <- diag(nrow=nrow(D_r))
  D_U <- diag(as.vector(Uhat))
  SigUstar <- (t(thetahatinv) %*% (D_A - (t(B) %*% D_p)) %*% D_r %*% thetahatinv) + (t(thetahatinv) %*% D_U %*% (D_r-Ident) %*% thetahatinv)
  ones <- rep(1,nrow(SigUstar))
  varNhat_marking_sum <- t(ones) %*% SigUstar %*% ones
  varNhat_marking <- diag(SigUstar)

  out$Nhat_strata1=as.vector(Nhatstar)
  out$se_Nhat_strata1=sqrt(as.vector(varNhat_marking))
  out$Nhat_strata2=as.vector(Nhat)
  out$se_Nhat_strata2=sqrt(as.vector(varNhat))
  out$Nhat=Nhat_sum
  out$se_Nhat=sqrt(as.vector(varNhat_sum))
  return(out)
}
#
# Mdata <- matrix(c(59,30,1,45,280,38,0,42,25),nrow=3,ncol=3,byrow=T)
# NDarroch(n1counts=c(484,1468,399),n2counts=c(847,6616,2489),stratamat=Mdata)
#
# Mdata <- matrix(5,nrow=4,ncol=3,byrow=T)
# NDarroch(n1counts=c(30,30,30,30),n2counts=c(30,40,50),stratamat=Mdata)
