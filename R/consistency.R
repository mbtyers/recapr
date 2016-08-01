#' Consistency Tests for the Abundance Estimator, Partial Stratification
#' @description Conducts three chi-squared tests for the consistency of the
#'   Petersen-type abundance estimator.  These tests provide explore evidence
#'   against the second traditional assumption of the Petersen mark-recapture
#'   experiment: that equal capture probabilities exist in either the first or
#'   second sampling event, or that complete mixing occurs between events.
#'
#'   Typically, if any of these test p-values is greater than the significance
#'   level, use of a Petersen-type estimator is considered justified.  If all
#'   three tests give p-values below the significance level and no movement
#'   occurs between strata (and strata are the same between events), a
#'   stratified estimator may be used. If all three tests give p-values below
#'   the significance level and some movement between strata occurs, a partially
#'   stratified (Darroch-type) estimator must be used, such as \link{NDarroch}.
#'
#'   This function assumes stratification in both sampling events, and in
#'   different ways (by time, area, etc.)  If stratification was the same in
#'   both events such that individuals could not move from one strata to another
#'   (such as by size or gender), use of \link{strattest} is recommended.
#' @param n1 A vector of the total sample sizes in the first event, by
#'   strata.  For example, setting \code{n1=c(20,30,40)} would imply 20
#'   individuals captured and marked in stratum 1, 30 in stratum 2, and 40 in
#'   stratum 3.
#' @param n2 A vector of the total sample sizes in the second event, by
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
#' @param ... Additional arguments for \link[stats]{chisq.test}
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @return A list of class \code{"recapr_consistencytest"} with the following components:
#' \itemize{
#' \item{\code{test1_tab}} {The contingency table used for the first test}
#' \item{\code{test1_Xsqd}} {The chi-squared test statistic in the first test}
#' \item{\code{test1_df}} {The associated degrees of freedom in the first test}
#' \item{\code{test1_pval}} {The p-value returned from the first test}
#' \item{\code{test2_tab}} {The contingency table used for the second test}
#' \item{\code{test2_Xsqd}} {The chi-squared test statistic in the second test}
#' \item{\code{test2_df}} {The associated degrees of freedom in the second test}
#' \item{\code{test2_pval}} {The p-value returned from the second test}
#' \item{\code{test3_tab}} {The contingency table used for the third test}
#' \item{\code{test3_Xsqd}} {The chi-squared test statistic in the third test}
#' \item{\code{test3_df}} {The associated degrees of freedom in the third test}
#' \item{\code{test3_pval}} {The p-value returned from the third test}
#' }
#' @author Matt Tyers
#' @seealso \link{strattest}, \link{NDarroch}
#' @importFrom stats chisq.test
#' @examples
#' consistencytest(n1=c(15,12,6), n2=c(12,9,10,8),
#'    m2strata1=c(1,1,1,1,1,2,2,2,3,3),
#'    m2strata2=c(1,1,3,3,4,1,2,4,1,3),
#'    simulate.p.value=TRUE)
#'
#' mat <- matrix(c(30,15,1,0,22,15), nrow=2, ncol=3, byrow=TRUE)
#' consistencytest(n1=c(284,199), n2=c(347,3616,1489), stratamat=mat)
#' @export
consistencytest <- function(n1, n2, m2strata1=NULL, m2strata2=NULL, stratamat=NULL,...) {
  if(!is.null(stratamat)) {
    if((nrow(stratamat) != length(n1)) | (ncol(stratamat) != length(n2))) stop("Dimension mismatch in stratamat - must have rows corresponding to n1 strata and columns corresponding to n1 strata")
  }
  if(is.null(stratamat)) {
    if(is.null(m2strata1) | is.null(m2strata2)) stop("recapture strata must be specified, either with stratamat or m2strata1 and m2strata2 together")
    if(any(!is.numeric(m2strata1))) stop("m2strata1 values must be positive integers")
    if(any(!is.numeric(m2strata2))) stop("m2strata2 values must be positive integers")
    if(sum(m2strata1-round(m2strata1))!=0) stop("m2strata1 values must be positive integers")
    if(sum(m2strata2-round(m2strata2))!=0) stop("m2strata2 values must be positive integers")
    if(any(m2strata1<1)) stop("m2strata1 values must be positive integers")
    if(any(m2strata2<1)) stop("m2strata2 values must be positive integers")
    if(max(m2strata1)>length(n1)) stop("strata values larger than count vector in event 1")
    if(max(m2strata2)>length(n2)) stop("strata values larger than count vector in event 2")
    if(length(m2strata1)!=length(m2strata2)) stop("m2strata1 and m2strata2 must be the same length (each element corresponding to an individual)")
    stratamat <- matrix(NA,nrow=length(n1),ncol=length(n2))
    for(i in 1:length(n1)) {
      for(j in 1:length(n2)) {
        stratamat[i,j] <- sum(m2strata1==i & m2strata2==j)
      }
    }
  }
  # stratamat <- table(m2strata1,m2strata2)
  mix_test_table <- cbind(stratamat,n1 - rowSums(stratamat))
  rownames(mix_test_table) <- 1:length(n1)
  colnames(mix_test_table) <- c(1:length(n2),"not recaptured")
  # cat('\n',"Mixing test",'\n')
  # cat("H0: Movement probabilities from stratum i to stratum j are the same among sections (all theta_ij = theta_j)",'\n','\n')
  # print(mix_test_table)
  suppressWarnings(x<-chisq.test(mix_test_table,...=...))
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  m22 <- colSums(stratamat)
  unmarked <- n2-m22
  eq_prop_table <- rbind(m22,unmarked)
  rownames(eq_prop_table) <- c("recaptured","unmarked")
  colnames(eq_prop_table) <- 1:length(n2)
  # cat('\n',"Equal proportions test",'\n')
  # cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  # print(eq_prop_table)
  suppressWarnings(y<-chisq.test(eq_prop_table,...=...))
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  m21 <- rowSums(stratamat)
  notrecap <- n1-m21
  complete_mix_table <- rbind(m21,notrecap)
  rownames(complete_mix_table) <- c("recaptured","not recaptured")
  colnames(complete_mix_table) <- 1:length(n1)
  # cat('\n',"Complete mixing test",'\n')
  # cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  # print(complete_mix_table)
  suppressWarnings(z<-chisq.test(complete_mix_table,...=...))
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  if(any(c(as.vector(x$expected),as.vector(y$expected),as.vector(z$expected))<5)) warning("Expected counts < 5, chi-squared approximation may be inaccurate.")

  out <- list(test1_tab=mix_test_table, test1_Xsqd=unname(x$statistic), test1_df=unname(x$parameter), test1_pval=unname(x$p.value),
              test2_tab=eq_prop_table, test2_Xsqd=unname(y$statistic), test2_df=unname(y$parameter), test2_pval=unname(y$p.value),
              test3_tab=complete_mix_table, test3_Xsqd=unname(z$statistic), test3_df=unname(z$parameter), test3_pval=unname(z$p.value))
  class(out) <- "recapr_consistencytest"
  return(out)
}


#' Print method for consistency test
#' @description Print method for consistency test
#' @param x Output from \code{consistencytest()}
#' @param ... additional print arguments
#' @author Matt Tyers
#' @method print recapr_consistencytest
#' @export
print.recapr_consistencytest <- function(x, ...) {
  cat("MIXING TEST",'\n')
  cat("H0: Movement probabilities from stratum i to stratum j are the same among sections (all theta_ij = theta_j)",'\n','\n')
  print(x$test1_tab)
  cat('\n',"X-squared: ",x$test1_Xsqd,"  df: ",x$test1_df,"  p-val: ",x$test1_pval,'\n','\n')

  cat("EQUAL PROPORTIONS TEST",'\n')
  cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  print(x$test2_tab)
  cat('\n',"X-squared: ",x$test2_Xsqd,"  df: ",x$test2_df,"  p-val: ",x$test2_pval,'\n','\n')

  cat("COMPLETE MIXING TEST",'\n')
  cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  print(x$test3_tab)
  cat('\n',"X-squared: ",x$test3_Xsqd,"  df: ",x$test3_df,"  p-val: ",x$test3_pval,'\n','\n')
  # print(x=...)
}

#consistencytest(n1=c(15,12,6), n2=c(12,9,10,8),m2strata1=c(1,1,1,1,1,2,2,2,3,3),m2strata2=c(1,1,3,3,4,1,2,4,1,3))

#' Power of Consistency Tests, Partial Stratification
#' @description Conducts power calculations of the chi-squared tests for the
#'   consistency of the Petersen-type abundance estimator, in a partial
#'   stratification setting, such as by time or geographic area.  In the case of
#'   partial stratification, individuals may move from one stratum to another
#'   between the first and second sampling events, and strata do not need to be
#'   the same between events.
#' @param n1 Vector of anticipated n1 counts (sample size in the first event),
#'   each element corresponding to one stratum.
#' @param n2 Vector of anticipated n2 counts (sample size in the second event),
#'   each element corresponding to one stratum.
#' @param pmat Matrix of assumed movement probabilities between strata, with
#'   rows corresponding to first-event strata and columns corresponding to
#'   second-event strata, and an additional column corresponding to the
#'   probability of NOT being recaptured in the second event.  Values will be
#'   standardized by row, that is, by first-event strata.  See note on usage
#'   below.
#' @param alpha Significance level for testing.  Defaults to \code{0.05}
#' @param sim Whether to conduct power calculation by simulation as well as
#'   Cohen's method.  Defaults to \code{TRUE}.
#' @param nsim Number of simulations if \code{sim} is \code{TRUE}.  Defaults to
#'   \code{10000}.
#' @return An object of class \code{"recapr_consistencypow"} with the following
#'   components: \itemize{ \item{\code{pwr1_c}} {Power of the first test,
#'   according to Cohen's method} \item{\code{pwr1_sim}} {Power of the first
#'   test, from simulation} \item{\code{ntest1}} {The sample size used for the
#'   first test} \item{\code{p0test1}} {The null-hypothesis probabilities for
#'   the first test} \item{\code{p1test1}} {The alt-hypothesis probabilities for
#'   the first test} \item{\code{pwr2_c}} {Power of the second test, according
#'   to Cohen's method} \item{\code{pwr2_sim}} {Power of the second test, from
#'   simulation} \item{\code{ntest2}} {The sample size used for the second test}
#'   \item{\code{p0test2}} {The null-hypothesis probabilities for the second
#'   test} \item{\code{p1test2}} {The alt-hypothesis probabilities for the
#'   second test} \item{\code{pwr3_c}} {Power of the third test, according to
#'   Cohen's method} \item{\code{pwr3_sim}} {Power of the third test, from
#'   simulation} \item{\code{ntest3}} {The sample size used for the third test}
#'   \item{\code{p0test3}} {The null-hypothesis probabilities for the third
#'   test} \item{\code{p1test3}} {The alt-hypothesis probabilities for the third
#'   test} \item{\code{alpha}} {The significance level used} }
#' @author Matt Tyers
#' @note The movement probability matrix specified in \code{pmat} is considered
#'   conditional on each row, that is, first-event strata, with columns
#'   corresponding to second-event strata and the final column specifying the
#'   probability of not being recaptured in the second event.  Values do not
#'   need to sum to one for each row, but will be standardized by the function
#'   to sum to one.
#'
#'   A \code{pmat} with a first row equal to \code{(0.05, 0.1, 0.15, 0.7)} would
#'   imply a 5 percent chance that individuals captured in the first-event
#'   strata 1 will be recaptured in second-event strata 1, and a 70 percent
#'   chance that individuals captured in the first-event strata 1 will not be
#'   recaptured in event 2.
#'
#'   Because of the row-wise scaling, specifying a row equal to \code{(0.05,
#'   0.1, 0.15, 0.7)} would be equivalent to that row having values \code{(1, 2, 3, 14)}.
#' @importFrom stats pchisq
#' @importFrom stats qchisq
#' @importFrom stats rmultinom
#' @references Cohen, J. (1988). Statistical power analysis for the behavioral
#'   sciences (2nd ed.). Hillsdale,NJ: Lawrence Erlbaum.
#'
#'   Code adapted from the 'pwr' package: Stephane Champely (2015). pwr: Basic
#'   Functions for Power Analysis. R package version 1.1-3.
#'   https://CRAN.R-project.org/package=pwr
#' @seealso \link{consistencytest}, \link{NDarroch}
#' @examples
#' mat <- matrix(c(4,3,2,1,10,3,4,3,2,10,2,3,4,3,10,1,2,3,4,10),
#'     nrow=4, ncol=5, byrow=TRUE)
#' powconsistencytest(n1=c(50,50,50,50), n2=c(50,50,50,50), pmat=mat)
#'
#' mat <- matrix(c(4,3,2,1,10,4,3,2,1,10,4,3,2,1,10,4,3,2,1,10),
#'     nrow=4, ncol=5, byrow=TRUE)
#' powconsistencytest(n1=c(50,50,50,50), n2=c(50,50,50,50), pmat=mat)
#'
#' mat <- matrix(c(1,1,1,1,10,2,2,2,2,10,3,3,3,3,10,4,4,4,4,10),
#'     nrow=4, ncol=5, byrow=TRUE)
#' powconsistencytest(n1=c(50,50,50,50), n2=c(50,50,50,50), pmat=mat)
#'
#' mat <- matrix(c(1,1,1,1,10,1,1,1,1,10,1,1,1,1,10,1,1,1,1,10),
#'     nrow=4, ncol=5, byrow=TRUE)
#' powconsistencytest(n1=c(50,50,50,50), n2=c(20,30,40,50), pmat=mat)
#'
#' mat <- matrix(c(1,1,1,1,5,1,1,1,1,8,1,1,1,1,10,1,1,1,1,15),
#'     nrow=4, ncol=5, byrow=TRUE)
#' powconsistencytest(n1=c(50,50,50,50), n2=c(50,50,50,50), pmat=mat)
#' @export
powconsistencytest <- function(n1,n2,pmat,alpha=0.05,sim=TRUE,nsim=10000) {
  if(length(n1) != nrow(pmat)) stop("pmat must have the same number of rows as first-event strata")
  if(length(n2) != (ncol(pmat)-1)) stop("pmat must have the same number of columns as second-event strata, plus one")

  nstrata1 <- length(n1)
  nstrata2 <- length(n2)

  mvtmat_rowcond <- pmat/matrix(rowSums(pmat),nrow=nrow(pmat),ncol=ncol(pmat),byrow=F)         # stratum 1->2 movement probabilities (conditional on row)
  mvtmat_expected <- mvtmat_rowcond*matrix(n1,nrow=nrow(pmat),ncol=ncol(pmat),byrow=F)         # stratum 1->2 expected counts
  probmat_expected <- rbind(mvtmat_expected,c((n2-colSums(mvtmat_expected[,1:length(n2)])),0)) # expected counts for all outcomes (except for not captured in either event)
  # probmat_expected[nrow(probmat_expected),ncol(probmat_expected)] <- N-sum(probmat_expected)
  # probmat_all <- probmat_expected/sum(probmat_expected)
  if(any(probmat_expected<0)) {
    stop("Negative expected counts generated")
    rownames(probmat_expected) <- 1:nstrata1
    colnames(probmat_expected) <- c(1:nstrata2,"not recap")
    print(probmat_expected)
  }

  #   p1 <- n1/N1
  #   p2 <- n2/N2
  #   pboth <- sum(p1)*sum(p2)
  #   probmat_both <- pmat/sum(pmat)*pboth
  #   probmat_all <- rbind(cbind(probmat_both,p1-rowSums(probmat_both)),c(p2-colSums(probmat_both),((1-sum(p1))*(1-sum(p2)))))
  #   if(any(p1<=rowSums(probmat_both))) stop("Negative first-event probabilities generated")
  #   if(any(p2<=colSums(probmat_both))) stop("Negative second-event probabilities generated")
  #
  #   probs_all <- as.vector(probmat_all)

  # probmat_mvt <- pmat/matrix(rowSums(pmat),nrow=nrow(pmat),ncol=ncol(pmat),byrow=F)                  # stratum 1->2 movement probabilities (conditional on row)
  # p1not <- 1-(sum(n1)/sum(N1))
  # p2not <- 1-(sum(n2)/sum(N2))
  # probmat_mvtall <- cbind(probmat_mvt*(1-p1not),rep(p1not,nstrata1))
  # expected_mvt <- probmat_mvtall*matrix(n1,nrow=nstrata1,ncol=(nstrata2+1),byrow=F)

  # pboth <- sum(p1)*sum(p2)
  # probmat_both <- pmat/sum(pmat)*pboth

  # probmat_mvt <- pmat/matrix(rowSums(pmat),nrow=nrow(pmat),ncol=ncol(pmat),byrow=F)                  # stratum 1->2 movement probabilities (conditional on row)
  # probmat_mvtall <- cbind((probmat_mvt*matrix(p1,nrow=nrow(pmat),ncol=ncol(pmat),byrow=F)),(1-p1))   # stratum 1->2 movement probabilities, including not recaptured (conditional on row)
  # probmat_both_m2 <- probmat_mvt*matrix(pN1*p1/sum(p1_crossclass),nrow=nrow(pmat),ncol=ncol(pmat),byrow=F)                    # total cross-classification probabilities (conditional on capture-recapture)
  # probmat_both <- probmat_both_m2*pboth                   # total cross-classification probabilities

  # probmat_all <- rbind(cbind(probmat_both,p1-rowSums(probmat_both)),c(p2-colSums(probmat_both),((1-sum(p1))*(1-sum(p2)))))
  # probmat_all <- rbind(cbind(probmat_both,p1_crossclass-rowSums(probmat_both)),c(p2_crossclass-colSums(probmat_both),((1-sum(p1_crossclass))*(1-sum(p2_crossclass)))))

  # if(any(p1<=rowSums(probmat_both))) stop("Negative first-event probabilities generated")
  # if(any(p2<=colSums(probmat_both))) stop("Negative second-event probabilities generated")

  # if(any(p1_crossclass<=rowSums(probmat_both))) stop("Negative first-event probabilities generated")
  # if(any(p2_crossclass<=colSums(probmat_both))) stop("Negative second-event probabilities generated")

  # probs_all <- as.vector(probmat_all)

  # mixing test - Cohen
  probmat_mixing <- probmat_expected[1:nstrata1,]/sum(probmat_expected[1:nstrata1,])
  probmat_mixingnull <- rowSums(probmat_mixing) %*% t(colSums(probmat_mixing))
  w1 <- sqrt(sum((probmat_mixing - probmat_mixingnull)^2/probmat_mixingnull))
  pwr1_c <- pchisq(q=qchisq(alpha,df=(nstrata1-1)*nstrata2,lower.tail=F),df=(nstrata1-1)*nstrata2,ncp=(sum(n2)*(w1^2)),lower.tail=F)

  # Equal proportions test - Cohen
  p1a <- rbind(colSums(probmat_expected[1:nstrata1,1:nstrata2]),probmat_expected[(nstrata1+1),1:nstrata2])
  p1a <- p1a/sum(p1a)
  p1anull <- rowSums(p1a) %*% t(colSums(p1a))
  w2 <- sqrt(sum((p1a - p1anull)^2/p1anull))
  pwr2_c <- pchisq(q=qchisq(alpha,df=(nstrata2-1),lower.tail=F),df=(nstrata2-1),ncp=(sum(n2)*(w2^2)),lower.tail=F)


  # complete mixing test - Cohen
  p2a <- rbind(rowSums(probmat_expected[1:nstrata1,1:nstrata2]),probmat_expected[1:nstrata1,(nstrata2+1)])
  p2a <- p2a/sum(p2a)
  p2anull <- rowSums(p2a) %*% t(colSums(p2a))
  w3 <- sqrt(sum((p2a - p2anull)^2/p2anull))
  pwr3_c <- pchisq(q=qchisq(alpha,df=(nstrata1-1),lower.tail=F),df=(nstrata1-1),ncp=(sum(n1)*(w3^2)),lower.tail=F)

  # simulation   ------------ figure out a more efficient way to do this!
  # reject1 <- reject2 <- reject3 <- rep(NA,nsim)
  # capmat <- probmat_mixing*NA
  # for(i in 1:nsim) {
  #   # mixing test
  #   test1mat <- matrix(rmultinom(1,sum(n1),as.vector(probmat_mixing)),nrow=nstrata1,ncol=nstrata2+1)
  #
  #   # equal proportions test
  #   test2mat <- matrix(rmultinom(1,sum(n2),as.vector(p1a)),nrow=2,ncol=nstrata2)
  #
  #   # complete mixing test
  #   test3mat <- matrix(rmultinom(1,sum(n1),as.vector(p2a)),nrow=2,ncol=nstrata1)
  #
  #   suppressWarnings(reject1[i] <- chisq.test(test1mat)$p.value<alpha)
  #   suppressWarnings(reject2[i] <- chisq.test(test2mat)$p.value<alpha)
  #   suppressWarnings(reject3[i] <- chisq.test(test3mat)$p.value<alpha)
  # }

  # pwr1_sim <- mean(reject1,na.rm=T)
  # pwr2_sim <- mean(reject2,na.rm=T)
  # pwr3_sim <- mean(reject3,na.rm=T)

  if(sim) {
    simmat1 <- rmultinom(n=nsim,size=sum(n1),prob=as.vector(probmat_mixing))
    ids <- expand.grid(1:nstrata1,1:(nstrata2+1))
    rowid <- unname(ids[,1])
    colid <- unname(ids[,2])
    expmat1 <- NA*simmat1
    for(i in 1:nstrata1) {
      for(j in 1:(nstrata2+1)) {
        expmat1[rowid==i & colid==j,] <- colSums(simmat1[rowid==i,])*colSums(simmat1[colid==j,])/sum(n1)
      }
    }
    X2_1 <- colSums(((simmat1-expmat1)^2)/expmat1)
    pwr1_sim2 <- mean(pchisq(X2_1,df=(nstrata1-1)*nstrata2,lower.tail=F)<alpha,na.rm=T)

    simmat2 <- rmultinom(n=nsim,size=sum(n2),prob=as.vector(p1a))
    ids <- expand.grid(1:2,1:nstrata2)
    rowid <- unname(ids[,1])
    colid <- unname(ids[,2])
    expmat2 <- NA*simmat2
    for(i in 1:2) {
      for(j in 1:nstrata2) {
        expmat2[rowid==i & colid==j,] <- colSums(simmat2[rowid==i,])*colSums(simmat2[colid==j,])/sum(n2)
      }
    }
    X2_2 <- colSums(((simmat2-expmat2)^2)/expmat2)
    pwr2_sim2 <- mean(pchisq(X2_2,df=(nstrata2-1),lower.tail=F)<alpha,na.rm=T)

    simmat3 <- rmultinom(n=nsim,size=sum(n1),prob=as.vector(p2a))
    ids <- expand.grid(1:2,1:nstrata1)
    rowid <- unname(ids[,1])
    colid <- unname(ids[,2])
    expmat3 <- NA*simmat3
    for(i in 1:2) {
      for(j in 1:nstrata1) {
        expmat3[rowid==i & colid==j,] <- colSums(simmat3[rowid==i,])*colSums(simmat3[colid==j,])/sum(n1)
      }
    }
    X2_3 <- colSums(((simmat3-expmat3)^2)/expmat3)
    pwr3_sim2 <- mean(pchisq(X2_3,df=(nstrata1-1),lower.tail=F)<alpha,na.rm=T)
  }

  rownames(probmat_mixing) <- rownames(probmat_mixingnull) <- 1:nstrata1
  colnames(probmat_mixing) <- colnames(probmat_mixingnull) <- c(1:nstrata2,"Not recap")

  rownames(p1a) <- rownames(p1anull) <- c("Recaptured","Unmarked")
  colnames(p1a) <- colnames(p1anull) <- 1:nstrata2

  rownames(p2a) <- rownames(p2anull) <- c("Marked","Unmarked")
  colnames(p2a) <- colnames(p2anull) <- 1:nstrata1

  if(sim) out <- list(pwr1_c=pwr1_c,pwr1_sim=pwr1_sim2,ntest1=sum(n1),p1test1=probmat_mixing,p0test1=probmat_mixingnull,pwr2_c=pwr2_c,pwr2_sim=pwr2_sim2,ntest2=sum(n2),p1test2=p1a,p0test2=p1anull,pwr3_c=pwr3_c,pwr3_sim=pwr3_sim2,ntest3=sum(n1),p1test3=p2a,p0test3=p2anull,alpha=alpha)
  else out <- list(pwr1_c=pwr1_c,ntest1=sum(n1),p1test1=probmat_mixing,p0test1=probmat_mixingnull,pwr2_c=pwr2_c,ntest2=sum(n2),p1test2=p1a,p0test2=p1anull,pwr3_c=pwr3_c,ntest3=sum(n1),p1test3=p2a,p0test3=p2anull,alpha=alpha)
  class(out) <- "recapr_consistencypow"
  return(out)
  # return(list(pwr1_c=pwr1_c,pwr1_sim=pwr1_sim,pwr1_sim2=pwr1_sim2,pwr2_c=pwr2_c,pwr2_sim=pwr2_sim,pwr2_sim2=pwr2_sim2,pwr3_c=pwr3_c,pwr3_sim=pwr3_sim,pwr3_sim2=pwr3_sim2))
  # return(list(pwr1_c=pwr1_c,pwr1_sim=pwr1_sim,pwr2_c=pwr2_c,pwr2_sim=pwr2_sim,pwr3_c=pwr3_c,pwr3_sim=pwr3_sim))
}
mat <- matrix(c(1,2,3,204,2,3,2,11,2,3,2,13),nrow=3,byrow=T)
powconsistencytest(pmat=mat,n1=c(12,20,18),n2=c(10,15,12),nsim=10000)


#' Print method for consistency test power
#' @description Print method for consistency test power
#' @param x Output from \code{powconsistencytest()}
#' @param ... additional print arguments
#' @author Matt Tyers
#' @method print recapr_consistencypow
#' @export
print.recapr_consistencypow <- function(x,...) {
  cat("MIXING TEST",'\n')
  cat("H0: Movement probabilities from stratum i to stratum j are the same among sections (all theta_ij = theta_j)",'\n')
  cat("Sample size (first event): ",x$ntest1,"   Significance level: ",x$alpha,'\n')
  cat('\n',"Null hypothesis cross-classification probabilities: ",'\n')
  print(round(x$p0test1,4))
  cat('\n',"Alternative hypothesis cross-classification probabilities: ",'\n')
  print(round(x$p1test1,4))
  cat('\n',"Null hypothesis expected counts: ",'\n')
  print(round(x$p0test1*x$ntest1,2))
  cat('\n',"Alternative hypothesis expected counts: ",'\n')
  print(round(x$p1test1*x$ntest1,2))
  cat('\n')
  if(!is.null(x$pwr1_sim)) cat("Power: ",x$pwr1_c,"   Power (from simulation): ",x$pwr1_sim,'\n','\n','\n')
  else cat("Power: ",x$pwr1_c,'\n','\n','\n')

  cat("EQUAL PROPORTIONS TEST",'\n')
  cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n')
  cat("Sample size (second event): ",x$ntest2,"   Significance level: ",x$alpha,'\n')
  cat('\n',"Null hypothesis capture probabilities: ",'\n')
  print(round(x$p0test2,4))
  cat('\n',"Alternative hypothesis capture probabilities: ",'\n')
  print(round(x$p1test2,4))
  cat('\n',"Null hypothesis expected counts: ",'\n')
  print(round(x$p0test2*x$ntest2,2))
  cat('\n',"Alternative hypothesis expected counts: ",'\n')
  print(round(x$p1test2*x$ntest2,2))
  cat('\n')
  if(!is.null(x$pwr2_sim)) cat("Power: ",x$pwr2_c,"   Power (from simulation): ",x$pwr2_sim,'\n','\n','\n')
  else cat("Power: ",x$pwr2_c,'\n','\n','\n')

  cat("COMPLETE MIXING TEST",'\n')
  cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n')
  cat("Sample size (first event): ",x$ntest3,"   Significance level: ",x$alpha,'\n')
  cat('\n',"Null hypothesis capture probabilities: ",'\n')
  print(round(x$p0test3,4))
  cat('\n',"Alternative hypothesis recapture probabilities: ",'\n')
  print(round(x$p1test3,4))
  cat('\n',"Null hypothesis expected counts: ",'\n')
  print(round(x$p0test3*x$ntest3,2))
  cat('\n',"Alternative hypothesis expected counts: ",'\n')
  print(round(x$p1test3*x$ntest3,2))
  cat('\n')
  if(!is.null(x$pwr3_sim)) cat("Power: ",x$pwr3_c,"   Power (from simulation): ",x$pwr3_sim,'\n','\n','\n')
  else cat("Power: ",x$pwr3_c,'\n','\n','\n')
}


#' Consistency Tests for the Abundance Estimator, Complete Stratification
#' @description Conducts two chi-squared tests for the consistency of the
#'   Petersen-type abundance estimator.  These tests provide explore evidence
#'   against equal capture probabilities in either the first or second sampling
#'   event.  Also conducts a third chi-squared test of unequal capture
#'   probabilities between sampling events for each stratum, in the case of
#'   small sample sizes (fewer than 100 in either sampling event and fewer than
#'   30 recaptures), which may be used to suggest unequal capture probabilities
#'   in either the first or second event.
#'
#'   Typically, if either of the first two test p-values is greater than the
#'   significance level, use of a Petersen-type estimator is considered
#'   justified.
#'
#'   If tests give evidence of unequal capture probabilities between strata, a
#'   stratified estimator should be used, such as \link{Nstrat}.
#'
#'   This function assumes stratification in both sampling events, such that
#'   individuals cannot move from one strata to another (such as by size or
#'   gender).  If movement between strata may occur (such as in the case of
#'   stratification by time or area), use of \link{consistencytest} is
#'   recommended.
#' @param n1 Vector of n1 counts (sample size in the first event), each element
#'   corresponding to one stratum.
#' @param n2 Vector of n2 counts (sample size in the second event), each element
#'   corresponding to one stratum.
#' @param m2 Vector of m2 counts (number of recaptures in the second event),
#'   each element corresponding to one stratum.
#' @param ... Additional arguments for \link[stats]{chisq.test}
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @return A list of class \code{"recapr_strattest"} with the following
#'   components: \itemize{ \item{\code{event1_table}} {The contingency table
#'   used for the first test} \item{\code{event1_Xsqd}} {The chi-squared test
#'   statistic in the first test} \item{\code{event1_df}} {The associated
#'   degrees of freedom in the first test} \item{\code{event1_pval}} {The
#'   p-value returned from the first test} \item{\code{event2_table}} {The
#'   contingency table used for the second test} \item{\code{event2_Xsqd}} {The
#'   chi-squared test statistic in the second test} \item{\code{event2_df}} {The
#'   associated degrees of freedom in the second test} \item{\code{event2_pval}}
#'   {The p-value returned from the second test} \item{\code{event1v2_table}}
#'   {The contingency table used for the third test} \item{\code{event1v2_Xsqd}}
#'   {The chi-squared test statistic in the third test}
#'   \item{\code{event1v2_df}} {The associated degrees of freedom in the third
#'   test} \item{\code{event1v2_pval}} {The p-value returned from the second
#'   third} }
#' @author Matt Tyers
#' @seealso \link{powstrattest}, \link{Nstrat}, \link{consistencytest}
#' @importFrom stats chisq.test
#' @examples
#' strattest(n1=c(100,100), n2=c(50,200), m2=c(20,15))
#' @export
strattest <- function(n1,n2,m2,...) {
  if((length(n1)!=length(n2)) | (length(n2)!=length(m2))) stop("input vectors must be equal length")
  unmarked <- n2-m2
  eq_prop_table <- rbind(m2,unmarked)
  rownames(eq_prop_table) <- c("recaptured","unmarked")
  colnames(eq_prop_table) <- 1:length(n2)
  # cat('\n',"Equal proportions test",'\n')
  # cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  # print(eq_prop_table)
  suppressWarnings(x<-chisq.test(eq_prop_table,...=...))
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  notrecap <- n1-m2
  complete_mix_table <- rbind(m2,notrecap)
  rownames(complete_mix_table) <- c("recaptured","not recaptured")
  colnames(complete_mix_table) <- 1:length(n1)
  # cat('\n',"Complete mixing test",'\n')
  # cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  # print(complete_mix_table)
  suppressWarnings(y<-chisq.test(complete_mix_table,...=...))
  # cat('\n',"X-squared: ",y$statistic,"  df: ",y$parameter,"  p-val: ",y$p.value,'\n','\n')

  event1v2_table <- rbind(n1,n2)
  rownames(event1v2_table) <- c("first event","second event")
  colnames(event1v2_table) <- 1:length(n1)
  suppressWarnings(z <- chisq.test(event1v2_table,...=...))

  if(any(c(as.vector(x$expected),as.vector(y$expected),as.vector(z$expected))<5)) warning("Expected counts < 5, chi-squared approximation may be inaccurate.")

  out<-list(event1_table=eq_prop_table, event1_Xsqd=unname(x$statistic), event1_df=unname(x$parameter), event1_pval=unname(x$p.value),
            event2_table=complete_mix_table, event2_Xsqd=unname(y$statistic), event2_df=unname(y$parameter), event2_pval=unname(y$p.value),
            event1v2_table=event1v2_table, event1v2_Xsqd=unname(z$statistic), event1v2_df=unname(z$parameter), event1v2_pval=unname(z$p.value))
  class(out) <- "recapr_strattest"
  return(out)
}


#' Print method for stratification test
#' @description Print method for stratification test
#' @param x Output from \code{strattest()}
#' @param ... additional print arguments
#' @author Matt Tyers
#' @method print recapr_strattest
#' @export
print.recapr_strattest <- function(x, ...) {
  cat('\n',"Equal proportions test",'\n')
  cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  print(x$event1_table)
  cat('\n',"X-squared: ",x$event1_Xsqd,"  df: ",x$event1_df,"  p-val: ",x$event1_pval,'\n','\n')

  cat('\n',"Complete mixing test",'\n')
  cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  print(x$event2_table)
  cat('\n',"X-squared: ",x$event2_Xsqd,"  df: ",x$event2_df,"  p-val: ",x$event2_pval,'\n','\n')

  cat('\n',"First vs. second event test",'\n')
  cat("H0: Probabilities of capture are the same between first and second events (all p_i equal to respective p_j)",'\n','\n')
  print(x$event1v2_table)
  cat('\n',"X-squared: ",x$event1v2_Xsqd,"  df: ",x$event1v2_df,"  p-val: ",x$event1v2_pval,'\n','\n')
  # print(x=...)
}

#strattest(n1=c(100,100),n2=c(20,200),m2=c(10,15))




#' Power of Consistency Tests, Complete Stratification
#' @description Conducts power calculations of the chi-squared tests for the
#'   consistency of the Petersen-type abundance estimator, in a complete
#'   stratification setting.
#' @param N Vector of total abundance, with each element corresponding to one
#'   stratum.
#' @param n1 Vector of anticipated n1 counts (sample size in the first event),
#'   each element corresponding to one stratum.
#' @param n2 Vector of anticipated n2 counts (sample size in the second event),
#'   each element corresponding to one stratum.
#' @param alpha Significance level for testing.  Defaults to \code{0.05}
#' @param sim Whether to conduct power calculation by simulation as well as
#'   Cohen's method.  Defaults to \code{TRUE}.
#' @param nsim Number of simulations if \code{sim} is \code{TRUE}.  Defaults to
#'   \code{100000}.
#' @return A list of three elements, each with class \code{"recapr_stratpow"}
#'   with the following components: \itemize{ \item{\code{prob}} {A vector of
#'   capture probabilities corresponding to the alternative hypothesis
#'   investigated} \item{\code{prob_null}} {A vector of capture probabilities
#'   corresponding to the null hypothesis (all probabilities equal)}
#'   \item{\code{n}} {The sample size used for the test} \item{\code{alpha}} {The significance level used for testing}
#'   \item{\code{power}} {The test power, calculated by Cohen's method}
#'   \item{\code{power_sim}} {The test power, calculated via simulation} }
#' @author Matt Tyers
#' @importFrom stats pchisq
#' @importFrom stats qchisq
#' @importFrom stats rbinom
#' @importFrom stats rhyper
#' @references Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale,NJ: Lawrence Erlbaum.
#'
#' Code adapted from the 'pwr' package:
#' Stephane Champely (2015). pwr: Basic Functions for Power Analysis. R
#' package version 1.1-3. https://CRAN.R-project.org/package=pwr
#' @seealso \link{strattest}, \link{Nstrat}
#' @examples
#' powstrattest(N=c(10000,20000), n1=c(1000,2000), n2=c(200,200))
#' @export
powstrattest <- function(N,n1,n2,alpha=0.05,sim=TRUE,nsim=100000) {
  if(length(N)!=length(n1) | length(N) != length(n2)) stop("N, n1, and n2 vectors must be the same length.")
  nstrata <- length(N)

  p1 <- n1/N
  p1null <- rep(sum(n1)/sum(N),nstrata)
  p2 <- n2/N
  p2null <- rep(sum(n2)/sum(N),nstrata)
  p3 <- n2/(n1+n2)
  p3null <- rep(sum(n2)/sum(n1+n2),nstrata)

  # ----------- Cohen method ------------ #
  # matrix of total probabilities
  p1a <- rbind((n1/sum(N)),((N-n1)/sum(N)))
  p2a <- rbind((n2/sum(N)),((N-n2)/sum(N)))
  p3a <- rbind(n1/sum(n1+n2),n2/sum(n1+n2))

  p1anull <- rowSums(p1a) %*% t(colSums(p1a))
  p2anull <- rowSums(p2a) %*% t(colSums(p2a))
  p3anull <- rowSums(p3a) %*% t(colSums(p3a))

  # effect sizes and power calculations
  w1 <- sqrt(sum((p1a - p1anull)^2/p1anull))
  pwr1_c <- pchisq(q=qchisq(alpha,df=(nstrata-1),lower.tail=F),df=(nstrata-1),ncp=(sum(n2)*(w1^2)),lower.tail=F)
  w2 <- sqrt(sum((p2a - p2anull)^2/p2anull))
  pwr2_c <- pchisq(q=qchisq(alpha,df=(nstrata-1),lower.tail=F),df=(nstrata-1),ncp=(sum(n1)*(w2^2)),lower.tail=F)
  w3 <- sqrt(sum((p3a - p3anull)^2/p3anull))
  pwr3_c <- pchisq(q=qchisq(alpha,df=(nstrata-1),lower.tail=F),df=(nstrata-1),ncp=(sum(n1+n2)*(w3^2)),lower.tail=F)

  # ------------ Simulation ------------- #
  n1sim <- n2sim <- m2sim <- matrix(NA,nrow=nsim,ncol=nstrata)
  for(j in 1:nstrata) {
    n1sim[,j] <- rbinom(n=nsim,size=N[j],prob=p1[j])
    n2sim[,j] <- rbinom(n=nsim,size=N[j],prob=p2[j])
    m2sim[,j] <- rhyper(nn=nsim,m=n1sim[,j],n=(N[j]-n1sim[,j]),k=n2sim[,j])
  }
  chi_p1 <- chi_p2 <- chi_p3 <- rep(NA,nsim)
  unmarked <- n2sim-m2sim
  notrecap <- n1sim-m2sim
  n1tot <- rowSums(n1sim)
  n2tot <- rowSums(n2sim)
  m2tot <- rowSums(m2sim)
  unmarkedtot <- rowSums(unmarked)
  notrecaptot <- rowSums(notrecap)

  m2expect1 <- m2tot*n2sim/n2tot
  unmarkedexpect <- unmarkedtot*n2sim/n2tot
  X2_1 <- rowSums((((m2sim-m2expect1)^2)/m2expect1)+(((unmarked-unmarkedexpect)^2)/unmarkedexpect))
  chi_p1 <- pchisq(q=X2_1,df=(nstrata-1),lower.tail=F)

  m2expect2 <- m2tot*n1sim/n1tot
  notrecapexpect <- notrecaptot*n1sim/n1tot
  X2_2 <- rowSums((((m2sim-m2expect2)^2)/m2expect2)+(((notrecap-notrecapexpect)^2)/notrecapexpect))
  chi_p2 <- pchisq(q=X2_2,df=(nstrata-1),lower.tail=F)

  nallsim <- n1sim+n2sim
  nalltot <- n1tot+n2tot
  n1expect <- nallsim*n1tot/nalltot
  n2expect <- nallsim*n2tot/nalltot
  X2_3 <- rowSums((((n1sim-n1expect)^2)/n1expect)+(((n2sim-n2expect)^2)/n2expect))
  chi_p3 <- pchisq(q=X2_3,df=(nstrata-1),lower.tail=F)

  pwr1_sim <- mean(chi_p1<=alpha,na.rm=T)
  pwr2_sim <- mean(chi_p2<=alpha,na.rm=T)
  pwr3_sim <- mean(chi_p3<=alpha,na.rm=T)

  # ------------- output ------------- #
  Test1 <- list(test="First-event capture probabilities",prob=p1,prob_null=p1null,n=sum(n2),alpha=alpha,power=pwr1_c)
  Test2 <- list(test="Second-event capture probabilities",prob=p2,prob_null=p2null,n=sum(n1),alpha=alpha,power=pwr2_c)
  Test3 <- list(test="First-event vs. second-event capture probabilities",prob=p3,prob_null=p3null,n=sum(n1+n2),alpha=alpha,power=pwr3_c)

  if(sim) {
    Test1$power_sim <- pwr1_sim
    Test2$power_sim <- pwr2_sim
    Test3$power_sim <- pwr3_sim
  }

  class(Test1) <- class(Test2) <- class(Test3) <- "recapr_stratpow"

  out <- list(Test1=Test1,Test2=Test2,Test3=Test3)
  return(out)
}


#' Print method for stratification test power
#' @description Print method for stratification test power
#' @param x Output from \code{powstrattest()}
#' @param ... additional print arguments
#' @author Matt Tyers
#' @method print recapr_stratpow
#' @export
print.recapr_stratpow <- function(x,...) {
  cat(x$test,'\n','\n')
  cat("Detecting capture probabilities for each strata: ",x$prob,'\n')
  cat("As compared to null values: ",x$prob_null,'\n')
  cat("sample size: ",x$n," and significance level: ",x$alpha,'\n','\n')
  cat("power: ",x$power,'\n')
  if(!is.null(x$power_sim)) cat("power (from simulation): ",x$power_sim,'\n')
}

# powstrattest(N=c(10000,20000),n1=c(1000,1000),n2=c(200,200))

#
# N <- 5000
# p <- 0.2
#
# # first, a vector of 10000 fish lengths
# length <- rnorm(N,mean=100,sd=10)
#
# # vectors of TRUE or FALSE depending on whether each fish is in M, C, or R
# # capture probability for M depends on length, capture probability for C does not
# inM <- rbinom(N, size=1, prob=p*(length-min(length))/(max(length)-min(length)))==1
# inC <- rbinom(N, size=1, prob=p)==1
# inR <- inM & inC
#
# # subsetting by whether each fish is in M, C, or R
# plot(density(length[inM]), col="red", main="",lwd=2)
# lines(density(length[inC]), col="blue",lwd=2)
# legend("topleft",legend=c("M","C"),lwd=2,col=c("red","blue"))
#
# ks.test(length[inC], length[inR])
# # p-value = 0.05035
#
# ks.test(length[inC & !inR], length[inR])
# # p-value = 0.02074


#
# nsim <- 10000000
# ntags <- 120
# nlocations <- 20 # each is 5% of population
# spawners <- rbinom(nsim,size=ntags,p=0.6)
# wheredidtheygo <- rmultinom(nsim,size=spawners,prob=rep(1,nlocations))
# probdetected <- rowMeans(wheredidtheygo>0)
# mean(probdetected)
