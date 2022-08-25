#' Calculate the stem age of a phylogeny
#'
#' @description Calculate the stem age of a rooted ultrametric phylogeny or set of phylogenies. This function can used only on stem crown phylogenie(s) and on binary or non binary phylogenie(s).
#'
#' @param phylo Object of class \code{phylo} or \code{multiPhylo}. A rooted ultrametric phylogeny of class \code{phylo} or a set of rooted ultrametric phylogenies of class \code{multiPhylo}. The rooted ultrametric phylogenie(s) can have polytomie(s) (i.e. non binary tree).
#'
#' @details This function will calculate the stem age of a phylogeny or a set of phylogenies.
#'
#' @return Returns the value of the stem age in the same unit as the branch length of the phylogenie(s) (typically in million of years).
#'
#' @author Sophia Lambert
#'
#' @export
#'
#' @seealso \code{\link{crown_age}}


stem_age <- function(phylo){
  if(!class(phylo)=="phylo"&!class(phylo)=="multiPhylo"){
    stop("the phylo object should be of class `phylo` or `multiPhylo`")
  }
  if(class(phylo)=="phylo"){
    phylo <- list(phylo)
  }
  if(any(sapply(seq_along(phylo), function(i)
    (!"root.edge" %in% names(phylo[[i]])|is.null(phylo[[i]][["root.edge"]]))))){
    stop("the phylo or multiPhylo object should have a stem age for each phylogenies at phylo$root.egde or phylo[[i]]$root.egde")
  }
  seqphy <- seq_along(phylo)
  len <- length(seqphy)
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
  stem.age <- sapply(seqphy, function(i)
    age[i]+phylo[[i]]$root.edge)
  return(stem.age)
}
