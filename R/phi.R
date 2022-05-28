#' Prior distribution on sampling probabilitie(s)
#'
#' @description Calculating the probability density function of the sampling probabilitie(s) \eqn{y} assuming a beta distribution or a uniform distribution from a to b.
#'
#' @param beta Logical. If TRUE, a beta distribution is assumed on the sampling probabilitie(s) \eqn{y}.
#' @param unif Logical. If TRUE, a uniform distribution from a to b is assumed on the sampling probabilitie(s) \eqn{y}.
#'
#' @details This function is a closure. This functions takes \code{beta} and \code{unif} as arguments and creates another function with arguments \code{x}, \code{a} and \code{b}.
#'
#' @return Returns an object of class \code{"function"}. This function takes the following arguments and will return the probability density of the sampling probabilitie(s) \eqn{y} assuming a beta distribution or a uniform distribution from a to b:
#' \describe{
#'   \item{x}{Numeric vector. The value of the sampling probabilitie(s) \eqn{y}.}
#'   \item{a}{Numeric. The value of \eqn{\alpha} or \eqn{a} respectively for the beta or the uniform distribution. This value cannot be negative for both distribution and cannnot exceed \eqn{b} for the uniform distribution.}
#'   \item{b}{Numeric. The value of \eqn{\beta} or \eqn{b} respectively for the beta or the uniform distribution. This value cannot be negative for both distribution, cannnot be inferior to \eqn{a} and cannot exceed 1 for the uniform distribution.}
#' }
#'
#' @author Sophia Lambert
#'
#' @export

#@seealso \code{\link{likelihood_bdRho}}, \code{\link{fitMCMC_bdRho}}

phi <- function(beta = F, unif = F){ # by default put a = 0 and b = 1
  if(beta == T & unif == T)
    stop("choose one unique prior distribution on the sampling fraction(s)") # check that should stop the code
  if(beta == F & unif == F)
    stop("choose one unique prior distribution on the sampling fraction(s)") # check that should stop the code
  if(beta == T){
    function(x, a, b){
      if(a < 0 | b < 0 )
        stop("a and b should be higher than 0") # check that should stop the code
      ifelse((a > 0) & (b > 0),
             x^(a-1)*(1-x)^(b-1)/
               (factorial(a-1)*factorial(b-1)/
                  factorial(a+b-1)), 0)
    }
  }
  else if(unif == T){
    function(x, a, b){
      if(a < 0 | b < 0 )
        stop("a and b should be higher than 0") # check that should stop the code
      ifelse(a <= x & x <= b & 0 <= a & b <= 1 & b > 0 & b > a, 1/(b-a), 0)
    }
  }
}
