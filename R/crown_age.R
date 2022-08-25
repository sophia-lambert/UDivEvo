#' Calculate the crown age of a phylogeny
#'
#' @description Calculate the crown age of a rooted ultrametric phylogeny or set of phylogenies. This function can used on stem or crown phylogenie(s) and on binary or non binary phylogenie(s).
#'
#' @param phylo Object of class \code{phylo} or \code{multiPhylo}. A rooted ultrametric phylogeny of class \code{phylo} or a set of rooted ultrametric phylogenies of class \code{multiPhylo}. The rooted ultrametric phylogenie(s) can have polytomie(s) (i.e. non binary tree).
#'
#' @details This function will calculate the crown age of a phylogeny or a set of phylogenies.
#'
#' @return Returns the value of the crown age in the same unit as the branch length of the phylogenie(s) (typically in million of years).
#'
#' @author Sophia Lambert
#'
#' @export
#'
#' @seealso \code{\link{stem_age}}


crown_age <- function(phylo){
  if(!class(phylo)=="phylo"&!class(phylo)=="multiPhylo"){
    stop("the phylo object should be of class `phylo` or `multiPhylo`")
  }
  if(class(phylo)=="phylo"){
    phylo <- list(phylo)
  }
  seqphy <- seq_along(phylo)
  len <- length(seqphy)
  ## put a warning message if seqphy is not the same length as tot_time, y, r and epsi, yl
  nbtips <- sapply(phylo, FUN = ape::Ntip)
  from_past <- lapply(seqphy, function (i)
    cbind(phylo[[i]]$edge, picante::node.age(phylo[[i]])$ages))
  ages <- lapply(seqphy, function(i)
    rbind(from_past[[i]][, 2:3], c(nbtips[i] + 1, 0)))
  ages <- lapply(seqphy, function(i)
    ages[[i]][order(ages[[i]][, 1]), ])
  counts <- lapply(seqphy, function (i)
    table(phylo[[i]]$edge[,1])) # handling polytomies
  polytom_nodes <- lapply(seqphy, function (i)
    as.numeric(names(counts[[i]])[counts[[i]] > 2]))
  polytomTimes <- lapply(seqphy, function (i)
    counts[[i]][counts[[i]] > 2]-1)
  ages <- lapply(seqphy, function (i)
    cbind(ages[[i]], 1))
  ages <- lapply(seqphy, function(i){
    ages[[i]][polytom_nodes[[i]],3] <- polytomTimes[[i]]; ages[[i]]})
  ages_polytom <- lapply(seqphy, function (i)
    as.data.frame(lapply(as.data.frame(ages[[i]]), rep, ages[[i]][,3])))
  age <- sapply(seqphy, function(i)
    max(ages_polytom[[i]][, 2]))
  return(age)
}
