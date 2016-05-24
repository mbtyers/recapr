#' Chapman Estimator
#' @description Calculates the value of the Chapman estimator for abundance in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.  The Chapman estimator (Chapman modification of the Petersen
#'   estimator) typically outperforms the Petersen estimator, even though the
#'   Peterson estimator is the MLE.
#' @param n1 Number of individuals captured and marked in the first sample.
#'   This may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The value of the Chapman estimator, calculated as (n1+1)*(n2+1)/(m2+1) - 1
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @references Chapman, D.G. (1951). Some properties of the hypergeometric distribution with applications to zoological censuses.  \emph{Univ. Calif. Public. Stat.} \strong{1}, 131-60.
#' @seealso \link{NPetersen}, \link{NBailey}, \link{vChapman}, \link{seChapman}, \link{rChapman}, \link{pChapman},
#'   \link{powChapman}, \link{ciChapman}
#' @examples
#' NChapman(n1=100, n2=100, m2=20)
#' @export
NChapman <- function(n1, n2, m2) (n1+1)*(n2+1)/(m2+1)-1



#' Petersen Estimator
#' @description Calculates the value of the Petersen estimator for abundance in
#'   a mark-recapture experiment, with given values of sample sizes and number
#'   of recaptures.  The Petersen estimator is the MLE, but is typically
#'   outperformed by the Chapman estimator.
#' @param n1 Number of individuals captured and marked in the first sample.
#'   This may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The value of the Petersen estimator, calculated as n1*n2/m2
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{NBailey}, \link{vPetersen}, \link{sePetersen}, \link{rPetersen}, \link{pPetersen},
#'   \link{powPetersen}, \link{ciPetersen}
#' @examples
#' NPetersen(n1=100, n2=100, m2=20)
#' @export
NPetersen <- function(n1, n2, m2) (n1)*(n2)/(m2)



#' Bailey Estimator
#' @description Calculates the value of the Bailey estimator for abundance in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.  The Bailey estimator assumes a binomial probability model in
#'   the second sampling event (i.e. sampling with replacement), rather than the
#'   hypergeometric model assumed by the Petersen and Chapman estimators.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The value of the Bailey estimator, calculated as n1*(n2+1)/(m2+1)
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @references Bailey, N.T.J. (1951). On estimating the size of mobile populations from capture-recapture data.  \emph{Biometrika} \strong{38}, 293-306.
#' @references Bailey, N.T.J. (1952). Improvements in the interpretation of recapture data.  \emph{J. Animal Ecol.} \strong{21}, 120-7.
#' @seealso \link{NPetersen}, \link{NChapman}, \link{vBailey}, \link{seBailey}, \link{rBailey}, \link{pBailey},
#'   \link{powBailey}, \link{ciBailey}
#' @examples
#' NBailey(n1=100, n2=100, m2=20)
#' @export
NBailey <- function(n1, n2, m2) (n1)*(n2+1)/(m2+1)



#' Estimated Variance of the Chapman Estimator
#' @description Calculates the estimated variance of the Chapman estimator in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The estimate variance of the Chapman estimator, calculated as
#'   (n1+1)*(n2+1)*(n1-m2)*(n2-m2)/((m2+2)*(m2+1)^2)
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{seChapman}, \link{rChapman}, \link{pChapman},
#'   \link{powChapman}, \link{ciChapman}
#' @examples
#' vChapman(n1=100, n2=100, m2=20)
#' @export
vChapman <- function(n1, n2, m2) (n1+1)*(n2+1)*(n1-m2)*(n2-m2)/(m2+1)/(m2+1)/(m2+2)



#' Estimated Variance of the Petersen Estimator
#' @description Calculates the estimated variance of the Petersen estimator in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The estimate variance of the Petersen estimator, calculated as
#'   (n1^2)*n2*(n2-m2)/(m2^3)
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NPetersen}, \link{sePetersen}, \link{rPetersen}, \link{pPetersen},
#'   \link{powPetersen}, \link{ciPetersen}
#' @examples
#' vPetersen(n1=100, n2=100, m2=20)
#' @export
vPetersen <- function(n1, n2, m2) {
#   Nhat <- (n1)*(n2)/(m2)
#   vhat <- m2*Nhat/(Nhat-1)*(1-(n1/Nhat))*(1-(n2/Nhat))
#   return(vhat)
  (n1^2)*n2*(n2-m2)/(m2^3)
}



#' Estimated Variance of the Bailey Estimator
#' @description Calculates the estimated variance of the Bailey estimator in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The estimate variance of the Bailey estimator, calculated as
#'   (n1^2)*(n2+1)*(n2-m2)/(m2+1)/(m2+1)/(m2+2)
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NBailey}, \link{seBailey}, \link{rBailey}, \link{pBailey},
#'   \link{powBailey}, \link{ciBailey}
#' @examples
#' vBailey(n1=100, n2=100, m2=20)
#' @export
vBailey <- function(n1, n2, m2) (n1^2)*(n2+1)*(n2-m2)/(m2+1)/(m2+1)/(m2+2)


#' Standard Error of the Chapman Estimator
#' @description Calculates the standard error of the Chapman estimator in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The estimate variance of the Chapman estimator, calculated as
#'   sqrt((n1+1)*(n2+1)*(n1-m2)*(n2-m2)/((m2+2)*(m2+1)^2))
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NChapman}, \link{vChapman}, \link{rChapman}, \link{pChapman},
#'   \link{powChapman}, \link{ciChapman}
#' @examples
#' seChapman(n1=100, n2=100, m2=20)
#' @export
seChapman <- function(n1, n2, m2) sqrt(vChapman(n1, n2, m2))


#' Standard Error of the Petersen Estimator
#' @description Calculates the standard error of the Petersen estimator in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The estimate variance of the Petersen estimator, calculated as
#'   sqrt((n1^2)*n2*(n2-m2)/(m2^3))
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NPetersen}, \link{vPetersen}, \link{rPetersen}, \link{pPetersen},
#'   \link{powPetersen}, \link{ciPetersen}
#' @examples
#' sePetersen(n1=100, n2=100, m2=20)
#' @export
sePetersen <- function(n1, n2, m2) sqrt(vPetersen(n1, n2, m2))


#' Standard Error of the Bailey Estimator
#' @description Calculates the standard error of the Bailey estimator in a
#'   mark-recapture experiment, with given values of sample sizes and number of
#'   recaptures.
#' @param n1 Number of individuals captured and marked in the first sample. This
#'   may be a single number or vector of values.
#' @param n2 Number of individuals captured in the second sample.  This may be a
#'   single number or vector of values.
#' @param m2 Number of marked individuals recaptured in the second sample.  This
#'   may be a single number or vector of values.
#' @return The estimate variance of the Bailey estimator, calculated as
#'   sqrt((n1^2)*(n2+1)*(n2-m2)/(m2+1)/(m2+1)/(m2+2))
#' @note Any Petersen-type estimator (such as this) depends on a set of
#'   assumptions: \itemize{ \item  The population is closed; that is, that there
#'   are no births, deaths, immigration, or emigration between sampling events
#'   \item All individuals have the same probability of capture in one of the
#'   two events, or complete mixing occurs between events \item Marking in the
#'   first event does not affect probability of recapture in the second event
#'   \item Individuals do not lose marks between events \item All marks will be
#'   reported in the second event }
#' @author Matt Tyers
#' @seealso \link{NBailey}, \link{vBailey}, \link{rBailey}, \link{pBailey},
#'   \link{powBailey}, \link{ciBailey}
#' @examples
#' seBailey(n1=100, n2=100, m2=20)
#' @export
seBailey <- function(n1, n2, m2) sqrt(vBailey(n1, n2, m2))






