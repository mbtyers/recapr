#' Consistency Tests for the Abundance Estimator
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
#' consistencytest(n1counts=c(15,12,6), n2counts=c(12,9,10,8),
#'    m2strata1=c(1,1,1,1,1,2,2,2,3,3),
#'    m2strata2=c(1,1,3,3,4,1,2,4,1,3),
#'    simulate.p.value=TRUE)
#'
#' mat <- matrix(c(30,15,1,0,22,15), nrow=2, ncol=3, byrow=TRUE)
#' consistencytest(n1counts=c(284,199), n2counts=c(347,3616,1489), stratamat=mat)
#' @export
consistencytest <- function(n1counts, n2counts, m2strata1=NULL, m2strata2=NULL, stratamat=NULL,...) {
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
  }
  # stratamat <- table(m2strata1,m2strata2)
  mix_test_table <- cbind(stratamat,n1counts - rowSums(stratamat))
  rownames(mix_test_table) <- 1:length(n1counts)
  colnames(mix_test_table) <- c(1:length(n2counts),"not recaptured")
  # cat('\n',"Mixing test",'\n')
  # cat("H0: Movement probabilities from stratum i to stratum j are the same among sections (all theta_ij = theta_j)",'\n','\n')
  # print(mix_test_table)
  x<-chisq.test(mix_test_table,...=...)
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  m22 <- colSums(stratamat)
  unmarked <- n2counts-m22
  eq_prop_table <- rbind(m22,unmarked)
  rownames(eq_prop_table) <- c("marked","unmarked")
  colnames(eq_prop_table) <- 1:length(n2counts)
  # cat('\n',"Equal proportions test",'\n')
  # cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  # print(eq_prop_table)
  y<-chisq.test(eq_prop_table,...=...)
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  m21 <- rowSums(stratamat)
  notrecap <- n1counts-m21
  complete_mix_table <- rbind(m21,notrecap)
  rownames(complete_mix_table) <- c("recaptured","not recaptured")
  colnames(complete_mix_table) <- 1:length(n1counts)
  # cat('\n',"Complete mixing test",'\n')
  # cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  # print(complete_mix_table)
  z<-chisq.test(complete_mix_table,...=...)
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  out <- list(test1_tab=mix_test_table, test1_Xsqd=unname(x$statistic), test1_df=unname(x$parameter), test1_pval=unname(x$p.value),
              test2_tab=eq_prop_table, test2_Xsqd=unname(y$statistic), test2_df=unname(y$parameter), test2_pval=unname(y$p.value),
              test3_tab=complete_mix_table, test3_Xsqd=unname(z$statistic), test3_df=unname(z$parameter), test3_pval=unname(z$p.value))
  class(out) <- "recapr_consistencytest"
  return(out)
}


#' @method print recapr_consistencytest
#' @export
print.recapr_consistencytest <- function(x, ...) {
  cat('\n',"Mixing test",'\n')
  cat("H0: Movement probabilities from stratum i to stratum j are the same among sections (all theta_ij = theta_j)",'\n','\n')
  print(x$test1_tab)
  cat('\n',"X-squared: ",x$test1_Xsqd,"  df: ",x$test1_df,"  p-val: ",x$test1_pval,'\n','\n')

  cat('\n',"Equal proportions test",'\n')
  cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  print(x$test2_tab)
  cat('\n',"X-squared: ",x$test2_Xsqd,"  df: ",x$test2_df,"  p-val: ",x$test2_pval,'\n','\n')

  cat('\n',"Complete mixing test",'\n')
  cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  print(x$test3_tab)
  cat('\n',"X-squared: ",x$test3_Xsqd,"  df: ",x$test3_df,"  p-val: ",x$test3_pval,'\n','\n')
  # print(x=...)
}

#consistencytest(n1counts=c(15,12,6), n2counts=c(12,9,10,8),m2strata1=c(1,1,1,1,1,2,2,2,3,3),m2strata2=c(1,1,3,3,4,1,2,4,1,3))


#' Consistency Tests for the Abundance Estimator
#' @description Conducts two chi-squared tests for the consistency of the
#'   Petersen-type abundance estimator.  These tests provide explore evidence
#'   against equal capture probabilities in either the first or second sampling
#'   event.
#'
#'   Typically, if either of these test p-values is greater than the
#'   significance level, use of a Petersen-type estimator is considered
#'   justified.  If both tests give evidence of unequal capture probabilities
#'   between strata, a stratified estimator should be used, such as \link{Nstrat}.
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
#' @return A list of class \code{"recapr_strattest"} with the following components:
#' \itemize{
#' \item{\code{cap_table}} {The contingency table used for the first test}
#' \item{\code{cap_Xsqd}} {The chi-squared test statistic in the first test}
#' \item{\code{cap_df}} {The associated degrees of freedom in the first test}
#' \item{\code{cap_pval}} {The p-value returned from the first test}
#' \item{\code{recap_table}} {The contingency table used for the second test}
#' \item{\code{recap_Xsqd}} {The chi-squared test statistic in the second test}
#' \item{\code{recap_df}} {The associated degrees of freedom in the second test}
#' \item{\code{recap_pval}} {The p-value returned from the second test}
#' }
#' @author Matt Tyers
#' @seealso \link{strattest}, \link{Nstrat}
#' @importFrom stats chisq.test
#' @examples
#' strattest(n1=c(100,100), n2=c(50,200), m2=c(20,15))
#' @export
strattest <- function(n1,n2,m2) {
  if((length(n1)!=length(n2)) | (length(n2)!=length(m2))) stop("input vectors must be equal length")
  unmarked <- n2-m2
  eq_prop_table <- rbind(m2,unmarked)
  rownames(eq_prop_table) <- c("marked","unmarked")
  colnames(eq_prop_table) <- 1:length(n2)
  # cat('\n',"Equal proportions test",'\n')
  # cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  # print(eq_prop_table)
  x<-chisq.test(eq_prop_table)
  # cat('\n',"X-squared: ",x$statistic,"  df: ",x$parameter,"  p-val: ",x$p.value,'\n','\n')

  notrecap <- n1-m2
  complete_mix_table <- rbind(m2,notrecap)
  rownames(complete_mix_table) <- c("recaptured","not recaptured")
  colnames(complete_mix_table) <- 1:length(n1)
  # cat('\n',"Complete mixing test",'\n')
  # cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  # print(complete_mix_table)
  y<-chisq.test(complete_mix_table)
  # cat('\n',"X-squared: ",y$statistic,"  df: ",y$parameter,"  p-val: ",y$p.value,'\n','\n')

  out<-list(cap_table=eq_prop_table, cap_Xsqd=unname(x$statistic), cap_df=unname(x$parameter), cap_pval=unname(x$p.value),
            recap_table=complete_mix_table, recap_Xsqd=unname(y$statistic), recap_df=unname(y$parameter), recap_pval=unname(y$p.value))
  class(out) <- "recapr_strattest"
  return(out)
}


#' @method print recapr_strattest
#' @export
print.recapr_strattest <- function(x, ...) {
  cat('\n',"Equal proportions test",'\n')
  cat("H0: Equal probability of capture among n1 strata (all p_i equal)",'\n','\n')
  print(x$cap_table)
  cat('\n',"X-squared: ",x$cap_Xsqd,"  df: ",x$cap_df,"  p-val: ",x$cap_pval,'\n','\n')

  cat('\n',"Complete mixing test",'\n')
  cat("H0: Equal probability of recapture among n2 strata (all p_j equal)",'\n','\n')
  print(x$recap_table)
  cat('\n',"X-squared: ",x$recap_Xsqd,"  df: ",x$recap_df,"  p-val: ",x$recap_pval,'\n','\n')
  # print(x=...)
}

#strattest(n1=c(100,100),n2=c(20,200),m2=c(10,15))
