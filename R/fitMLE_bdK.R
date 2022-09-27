#' Maximum likelihood fit of the birth-death k model on a phylogeny
#'
#' @description Fits the birth-death k model (a constant-time homogeneous birth-death model under k-sampling) to a rooted ultrametric phylogeny using Maximum likelihood estimation. Precisely, it uses the Nelder-Mead algorithm to optimise the likelihood. The birth-death process is conditioned on the starting time of the process \code{tot_time} and the survival of the process at present time as well as having \eqn{k} extant sampled tips at present. This function can fit the birth-death k model on a stem or crown phylogeny. This function is specifically adapted for diversification analysis on phylogenies on which the sampling probability is unknown.
#'
#' @param phylo Object of class \code{phylo}. A rooted ultrametric phylogeny of class \code{phylo}. The rooted ultrametric phylogeny can have polytomie(s) (i.e. non binary tree).
#' @param tot_time Numeric. The stem or crown age (also called MRCA) of the phylogeny depending on the conditioning of the process specified (see \code{cond} argument) and the phylogeny used accordingly. The stem age of the phylogeny can be computed using max(TreeSim::getx(phylo))+phylo$root.edge (note that the phylo$root.edge needs to be known) and the crown age of the phylogeny can be computed using max(TreeSim::getx(phylo)).
#' @param r Numeric vector. The net diversification rate \eqn{r} starting value(s).
#' @param epsi Numeric vector. The turnover rate \eqn{\epsilon} starting value(s).
#' @param cond Character. Specifying the conditioning of the birth-death process. Two conditioning are available, either \code{cond = "crown"} (the default option) if the phylogeny used is a crown phylogeny or \code{cond = "stem"} if the phylogeny used is a stem phylogeny.
#' @param YULE Logical. If \code{TRUE}, the extinction rate \eqn{\mu} thus the turnover rate \eqn{\epsilon} are fixed to \eqn{0} and the net diversification rate \eqn{r} equals the speciation rate \eqn{\lambda}. If \code{FALSE} (the default option), the turnover rate \eqn{\epsilon} is not fixed to \eqn{0} and is thus inferred.
#' @param dt Numeric. If \code{dt = 0}, the integral on the sampling probability is computed using the R \code{stats::integrate} function. If \eqn{dt\ge0}, the integral of the sampling probability is performed manually using a piece-wise constant approximation. \code{dt} represents the length of the interval on which the function integrated is assumed to be constant. For manual integral, advised value of \code{dt} are \eqn{1e-3} to \eqn{1e-5}.
#' @param rel.tol Numeric. This argument is only used if \code{dt = 0}. This represents the relative accuracy requested when the integral is performed using the \code{stats::integrate} function. Typically \code{.Machine$double.eps^0.25} is used but a value of \eqn{1e-10} (the default value) has been tested and performs well.
#' @param tuned_dichotomy Logical. If \code{TRUE}, when the log likelihood of the model is equal to non finite value due to approximations, a dichotomy search is performed to find a tuning parameter that will be used for getting a finite value of the log likelihood. If \code{TRUE}, the log likelihood will take longer to calculate. Else if \code{FALSE}, no dichotomy search is performed; if the log likelihood is equal to non finite value due to approximations, the log likelihood will take this non finite value for the corresponding parameters.
#' @param brk Numeric. This argument is only used if \code{tuned_dichotomy = TRUE}. The number of steps used in the dichotomy search. Typically the value \eqn{200} is sufficient to avoid non finite values. In some case if the log likelihood is still equal to non finite value, the \code{brk} value \eqn{2000} will be required for more tuning but it will rarely take a larger value.
#'
#' @details This function will fit the birth-death k-sampling model and the function will infer the net diversification rate \eqn{r} and the turnover rate \eqn{\epsilon}. Note that the starting values should be adapted to the number of parameters (here \eqn{2}, \eqn{r} and \eqn{\epsilon}). This function is specifically intended to be used on phylogenies with unknown or highly uncertain global diversity estimates (the sampling probability is not known with accuracy). Note that the sampling probability is never estimated and that this function is not able to evaluate negative rates.
#'
#' @return Returns an object of class \code{MLE_bdK}. This \code{MLE_bdK} object is a list containing the name of the birth-death model performed, the log likelihood of the data knowing the parameters, the akaike information criterion corrected and the inferred parameters.
#'
#' @author Sophia Lambert
#'
#' @example /inst/examples/ExfitMLE_bdK.R
#'
#' @export
#'
#' @seealso \code{\link{likelihood_bdK}} and \code{\link{fitMCMC_bdK}}


fitMLE_bdK <- function(phylo, tot_time,
                        r, epsi,
                        cond = "crown",
                        YULE = FALSE, dt = 0,
                        rel.tol = 1e-10,
                        tuned_dichotomy = TRUE,
                        brk = 2000){

  # TO DO : rel.tol, if the inference stop by itself it means that the rel.tol is probably too low maybe add an option if the integrate is NULL and does not exist
  # same thing for save_inter with the checkpointing

  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"") # check stops

  nbtip <- ape::Ntip(phylo)
  from_past <- cbind(phylo$edge, picante::node.age(phylo)$ages)
  ages <- rbind(from_past[, 2:3], c(nbtip + 1, 0))
  counts <- table(phylo$edge[,1])
  polytom_nodes <- as.numeric(names(counts)[counts > 2])
  polytomTimes <- counts[counts > 2]-1
  ages <- cbind(ages, 1)
  polytom_nodes_position <- which(ages[,1] %in% polytom_nodes) # essential for the good position (because not ordered)
  ages[polytom_nodes_position,3] <- polytomTimes
  ages_polytom <- as.matrix(as.data.frame(lapply(as.data.frame(ages), rep, ages[,3]))[,1:2])
  colnames(ages_polytom)<-NULL
  age <- max(ages_polytom[,2])
  tj <- list()
  if (cond=="crown")
  {
    root <- 2
    ntot <- nbtip-1+nbtip
    from_past_polytom <- as.matrix(as.data.frame(lapply(as.data.frame(from_past), rep, ages[-length(ages[,3]),3])))
    select <- which(from_past_polytom[,1]==nbtip+1)[2]
    nbtips <- c()
    nbtips[1] <- length(which(ages_polytom[1:select-1,1]<(nbtip+1)))
    nbtips[2] <- length(which(ages_polytom[select:(ntot-1),1]<(nbtip+1)))
    out <-  which(ages_polytom[,1]%in%c(1:nbtip))
    ages_nodes <- ages_polytom[-out,]
    node <- list()
    node[[1]] <- 1:(nbtips[1]-1)
    node[[2]] <- nbtips[1]:(nbtips[1]+nbtips[2]-2)
    node[sapply(node, is.unsorted)] <- 0 # essential : if there is only one tip in a subphylo there is no tj for it
    tj[[1]] <- age - ages_nodes[node[[1]], 2]
    tj[[2]] <- age - ages_nodes[node[[2]], 2]
  }

  else if (cond=="stem")
  {
    root <- 1
    nbtips <- nbtip
    j <- 1:(as.numeric(nbtips) - 1)
    node <- as.numeric(nbtips) + j
    ages_polytom <- ages_polytom[order(ages_polytom[, 1]), ]
    tj[[1]] <- age - ages_polytom[node, 2]
  }

  likeli <- likelihood_bdK(tottime = tot_time, nbtips, tj,
                           root, dt = dt, rel.tol = rel.tol,
                           tuned_dichotomy = tuned_dichotomy,
                           brk = brk)

  # if we want to use the model with K sampling
  tun.init <- rep(log(1), root)
  if (YULE==FALSE)
  {
    if(is.null(r)|is.null(epsi))
      stop("\"r\" and \"epsi\" starting value should be specified") # checking if the stop actually stops

    startV <- c(r, epsi)
    p <- length(startV)

    optimLH <- function(init)
    {
      r <- init[1:length(p/2)]
      epsi <- init[(1 + length(p/2)):length(init)]
      tuningLH <- likeli(tun.init, seqphy.init = 1:root)
      bfLH <- tuningLH(div = r, turn = epsi)
      LH <- bfLH$logLik
      # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
      return(-LH)
    }
    names = c(paste("div rate", 1:length(p/2)), paste("epsi rate", 1:length(p/2)))
  }

  # if we want to use the model with K sampling and we force the turnover rate to be 0

  else # YULE = TRUE means fix.mu = TRUE and mu = 0 the d parameter becomes the lambda
  {
    message('extinction rate(s) and turnover rate(s) "epsi" are fixed to 0')

    if(is.null(r))
      stop("\"r\" starting value (becoming \"lambda\" when YULE = TRUE) should be specified") # checking if the stop actually stops

    startV <- c(r)
    p <- length(startV)

    optimLH <- function(init)
    {
      tuningLH <- likeli(tun.init, seqphy.init = 1:root)
      bfLH <- tuningLH(div = init, turn = 0)
      LH <- bfLH$logLik
      # if(tun_mem == TRUE) tun.init <<- bfLH$tun # IMPORTANT tun_mem does not exist anymore
      return(-LH)
    }
    names = c(paste("div rate", 1:length(p)))
  }

  temp <- stats::optim(startV, optimLH, method = "Nelder-Mead")
  df_params <- data.frame(temp$par, row.names = names)
  res <- list(model = "bd.K", LH = -temp$value,
              aicc=2 * temp$value + 2 * p + (2 * p * (p + 1)) / (nbtip - p - 1),
              parameters = df_params)

  class(res) <- "MLE_bdK"
  return(res)
}

