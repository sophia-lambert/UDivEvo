% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_bdK.R
\name{likelihood_bdK}
\alias{likelihood_bdK}
\title{Likelihood birth-death k function}
\usage{
likelihood_bdK(
  tottime,
  nbtips,
  tj,
  root,
  dt = dt,
  rel.tol = rel.tol,
  tuned_dichotomy = tuned_dichotomy,
  brk = brk
)
}
\arguments{
\item{tottime}{Numeric vector. The stem or crown age (also called MRCA) of the phylogeny depending on the conditioning of the process specified (see \code{root} argument) and the phylogeny used. The stem age of the phylogeny can be computed using max(TreeSim::getx(phylo))+phylo$root.edge (note that the phylo$root.edge needs to be known) and the crown age of the phylogeny can be computed using max(TreeSim::getx(phylo)).}

\item{nbtips}{Integer.  The number of extant sampled tips in the phylogeny typically noted \code{k}.}

\item{tj}{List of numeric vector. A list of atomic numeric vector where the vector is specifying the node depths of the phylogeny. The list is of length 1. The node depths of a phylogeny can be computed using \code{TreeSim::getx(phylo)[!TreeSim::getx(phylo)==tottime]} in a list. This code works both for stem and crown phylogenies. Note that the \code{tottime} is not contained in node depths here since it is used for conditioning the process. Since there are \eqn{k-1} internal nodes in a phylogeny and that the MRCA is not included in the node depths when conditioning in crown, the numeric vector is of length \eqn{k-2} if the phylogeny is in crown and of length \eqn{k-1} if it is in stem with \code{k} being the number of extant sampled tips in the phylogeny.}

\item{root}{Integer. Specifying the conditioning of the birth-death process. \code{root} takes the value \eqn{0} if the phylogeny used is a stem phylogeny to condition the process on the stem age. It takes the value \eqn{1} if the phylogeny is a crown phylogeny to condition the process on the crown age.}

\item{dt}{Numeric. If \code{dt = 0}, the integral on the sampling probability is computed using the R \code{stats::integrate} function. If \eqn{dt\ge0}, the integral of the sampling probability is performed manually using a piece-wise constant approximation. \code{dt} represents the length of the interval on which the function integrated is assumed to be constant. For manual integral, advised value of \code{dt} are \eqn{1e-3} to \eqn{1e-5}.}

\item{rel.tol}{Numeric. This represents the relative accuracy requested when the integral is performed using \code{stats::integrate} function. Typically, \code{.Machine$double.eps^0.25} is used but a value of \eqn{1e-10} has been tested and performs well.}

\item{tuned_dichotomy}{Logical. If \code{TRUE}, when the log likelihood is equal to non finite value due to approximations, a dichotomy search is performed to find a tuning parameter that will be used for getting a finite value of the log likelihood. If \code{TRUE}, the log likelihood will take longer to calculate. Else if \code{FALSE}, no dichotomy search is performed; if the log likelihood is equal to non finite value due to approximations, the function will return this non finite value.}

\item{brk}{Numeric. The number of steps used in the dichotomy search. Typically the value \eqn{200} is sufficient to avoid non finite values. In some case if the log likelihood is still equal to non finite value, the \code{brk} value \eqn{2000} will be required for more tuning but it will rarely take a larger value.}
}
\value{
Returns an object of class \code{function}. This function can take the following arguments (\code{tun.init} and \code{seqphy.init}) and will return another function. This last function will take the diversifications rate as arguments ((\code{div} and \code{turn})) and will return the value of the log likelihood in a list together with the tuning parameter value used:
\describe{
  \item{tun.init}{Numeric The initial tuning parameter value. Typically, it will take the value \eqn{log(1)}.}
  \item{seqphy.init}{Integer. Here \code{seqphy.init = 1} since it is the likelihood of a single phylogeny that is computed.}
  \item{div}{Numeric. The net diversification rate also called \eqn{r}.}
  \item{turn}{Numeric. The turnover rate also called \eqn{\epsilon}.}
}

Note that this functioon can be used on a set of phylogenies. See \code{\link{likelihood_bdRho}} to check how to adapt the arguments for a set of phylogenies.
}
\description{
Computes the log likelihood of a rooted ultrametric phylogeny under a constant-time homogeneous birth-death model and k-sampling. The birth-death process is conditioned on the starting time of the process \code{tottime} and the survival of the process at present time as well as having \eqn{k} extant sampled tips at present. This function can computes the log likelihood on a stem or crown phylogeny. This function is specifically adapted for diversification analysis on phylogenies on which the sampling probability is unknown.
}
\details{
This function is a closure, it takes all of the above as arguments and creates another function using the arguments described in the \code{Value} section. Note that when the phylogeny is in stem, the conditioning is done on the stem age and \eqn{k} number of extant tips while when the phylogeny is in crown, the conditioning is done on the crown age, \eqn{k_1} and \eqn{k_2} the number of extant tips of the sister clades diverging at the MRCA. This function is specifically intended to be used on phylogenies with unknown or highly uncertain global diversity estimates (the sampling probability is not known with accuracy). Note that the sampling probability is never estimated and that this function is not able to evaluate negative rates.
}
\seealso{
\code{\link{fitMCMC_bdK}} and \code{\link{likelihood_bdRho}}
}
\author{
Sophia Lambert
}
