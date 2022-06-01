#' Likelihood birth-death k function
#'
#' @description Computes the log likelihood of a rooted ultrametric phylogeny under a constant-time homogeneous birth-death model and k-sampling. The birth-death process is conditioned on the starting time of the process \code{tottime} and the survival of the process at present time as well as having \eqn{k} extant sampled tips at present. This function can computes the log likelihood on a stem or crown phylogeny. This function is specifically adapted for diversification analysis on phylogenies on which the sampling probability is unknown.
#'
#' @param tottime Numeric vector. The stem or crown age (also called MRCA) of the phylogeny depending on the conditioning of the process specified (see \code{root} argument) and the phylogeny used. The stem age of the phylogeny can be computed using max(TreeSim::getx(phylo))+phylo$root.edge (note that the phylo$root.edge needs to be known) and the crown age of the phylogeny can be computed using max(TreeSim::getx(phylo)).
#' @param nbtips Integer.  The number of extant sampled tips in the phylogeny typically noted \code{k}.
#' @param tj List of numeric vector. A list of atomic numeric vector where the vector is specifying the node depths of the phylogeny. The list is of length 1. The node depths of a phylogeny can be computed using \code{TreeSim::getx(phylo)[!TreeSim::getx(phylo)==tottime]} in a list. This code works both for stem and crown phylogenies. Note that the \code{tottime} is not contained in node depths here since it is used for conditioning the process. Since there are \eqn{k-1} internal nodes in a phylogeny and that the MRCA is not included in the node depths when conditioning in crown, the numeric vector is of length \eqn{k-2} if the phylogeny is in crown and of length \eqn{k-1} if it is in stem with \code{k} being the number of extant sampled tips in the phylogeny.
#' @param root Integer. Specifying the conditioning of the birth-death process. \code{root} takes the value \eqn{0} if the phylogeny used is a stem phylogeny to condition the process on the stem age. It takes the value \eqn{1} if the phylogeny is a crown phylogeny to condition the process on the crown age.
#' @param dt Numeric. If \code{dt = 0}, the integral on the sampling probability is computed using the R \code{stats::integrate} function. If \eqn{dt\ge0}, the integral of the sampling probability is performed manually using a piece-wise constant approximation. \code{dt} represents the length of the interval on which the function integrated is assumed to be constant. For manual integral, advised value of \code{dt} are \eqn{1e-3} to \eqn{1e-5}.
#' @param rel.tol Numeric. This represents the relative accuracy requested when the integral is performed using \code{stats::integrate} function. Typically, \code{.Machine$double.eps^0.25} is used but a value of \eqn{1e-10} has been tested and performs well.
#' @param tuned_dichotomy Logical. If \code{TRUE}, when the log likelihood is equal to non finite value due to approximations, a dichotomy search is performed to find a tuning parameter that will be used for getting a finite value of the log likelihood. If \code{TRUE}, the log likelihood will take longer to calculate. Else if \code{FALSE}, no dichotomy search is performed; if the log likelihood is equal to non finite value due to approximations, the function will return this non finite value.
#' @param brk Numeric. The number of steps used in the dichotomy search. Typically the value \eqn{200} is sufficient to avoid non finite values. In some case if the log likelihood is still equal to non finite value, the \code{brk} value \eqn{2000} will be required for more tuning but it will rarely take a larger value.
#'
#' @details This function is a closure, it takes all of the above as arguments and creates another function using the arguments described in the \code{Value} section. Note that when the phylogeny is in stem, the conditioning is done on the stem age and \eqn{k} number of extant tips while when the phylogeny is in crown, the conditioning is done on the crown age, \eqn{k_1} and \eqn{k_2} the number of extant tips of the sister clades diverging at the MRCA. This function is specifically intended to be used on phylogenies with unknown or highly uncertain global diversity estimates (the sampling probability is not known with accuracy). Note that the sampling probability is never estimated and that this function is not able to evaluate negative rates.
#'
#' @return Returns an object of class \code{function}. This function can take the following arguments (\code{tun.init} and \code{seqphy.init}) and will return another function. This last function will take the diversifications rate as arguments ((\code{div} and \code{turn})) and will return the value of the log likelihood in a list together with the tuning parameter value used:
#' \describe{
#'   \item{tun.init}{Numeric The initial tuning parameter value. Typically, it will take the value \eqn{log(1)}.}
#'   \item{seqphy.init}{Integer. Here \code{seqphy.init = 1} since it is the likelihood of a single phylogeny that is computed.}
#'   \item{div}{Numeric. The net diversification rate also called \eqn{r}.}
#'   \item{turn}{Numeric. The turnover rate also called \eqn{\epsilon}.}
#' }
#' @return Note that this functioon can be used on a set of phylogenies. See \code{\link{likelihood_bdRho}} to check how to adapt the arguments for a set of phylogenies.
#'
#' @author Sophia Lambert
#'
#' @export
#'
#' @seealso \code{\link{fitMCMC_bdK}} and \code{\link{likelihood_bdRho}}


likelihood_bdK <- function(tottime, nbtips, tj,
                           root, dt = dt, rel.tol = rel.tol,
                           tuned_dichotomy = tuned_dichotomy,
                           brk = brk){

  nbtip <- sum(nbtips)

  integr_Rho <- int_Rho(dt = dt)

  function(tun.init, seqphy.init){
    function(div, turn){

      tun = tun.init
      seqphy = seqphy.init
      if(any(div < 0) | any(turn < 0)){rv_fin = -Inf}
      else{
        like_tuning <- function(tun, seqphy){
          integr_loglik <- sapply(seqphy, function(i){
            log(integr_Rho(function(yj) sapply(yj, function (yj){
              exp(log(1-Pa(t = tottime, d = div, epsi = turn)*(1-yj))*-2 +
                    sum(log(Pnd(tj[[i]], d = div, epsi = turn, rho = yj)) +
                          tj[[i]]*div-log(div)) - tun[i]*(nbtips[i]-1))
            }), 0, 1, rel.tol = rel.tol))
          })
          return(integr_loglik)
        }

        rv = like_tuning(tun = tun, seqphy = seqphy)

        if(tuned_dichotomy == TRUE){
          try_int = 0
          tun <- tun[sapply(tj, length)>0]
          tun1 <- 10^308/(nbtips[sapply(tj, length)>0]-1) # will be the -Inf bound and to avoid having NaN (because we are trying to do -Inf+Inf)
          tun2 <- -10^308/(nbtips[sapply(tj, length)>0]-1) # will be the Inf bound and to avoid having NaN because we are trying to do Inf-Inf
          while (any(is.infinite(rv[sapply(tj, length)>0]))){
            rin  = rv[sapply(tj, length)>0] == Inf
            rmin = rv[sapply(tj, length)>0] == -Inf
            tun1[rmin] <- tun[rmin]
            tun2[rin] <- tun[rin]
            tun <- (tun1+tun2)/2
            if(any(round(tun1, digits = 10) == round(tun2, digits = 10))){ # describe that in the md
              break
              warning("tuning parameters too close to each other, the likelihood is hard to calculate")}

            wtr = rin | rmin
            seqphy <- which(wtr==TRUE)
            rv[sapply(tj, length)>0][wtr] = like_tuning(tun, which(sapply(tj, length)>0)[seqphy])

            try_int = try_int + 1
            if(try_int > brk) {
              break
              warning('could not calculate the likelihood with the bisection method up to ', brk,' splits \n')}
          }
        }

        else(rv)

        if(any(is.na(rv))){ # check if still useful
          print("NA detected, once passed in log this equals NaN")
          rv_fin <- -Inf
        }
        else if(is.na(sum(rv))){ # probably still useful
          print("NA detected probably because at least one of the two subtrees log likelihood phylo equals Inf and the other one equals -Inf thus the sum is NaN")
          rv_fin <- -Inf
        }
        else(rv_fin <- sum(rv)-
               sum(unlist(tj)*div-log(div)) + sum(tun*(nbtips-1)) +
               sum(log(nbtips)) +
               log(P1_a(t = tottime, d = div, epsi = turn))*root -
               log(Pa(t = tottime, d = div, epsi = turn)) * (nbtip - root))
      }
      resLH <- list(logLik = rv_fin, tuning = tun)
      return(resLH)
    }
  }
}
