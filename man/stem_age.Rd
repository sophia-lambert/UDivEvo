% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stem_age.R
\name{stem_age}
\alias{stem_age}
\title{Calculate the stem age of a phylogeny}
\usage{
stem_age(phylo)
}
\arguments{
\item{phylo}{Object of class \code{phylo} or \code{multiPhylo}. A rooted ultrametric phylogeny of class \code{phylo} or a set of rooted ultrametric phylogenies of class \code{multiPhylo}. The rooted ultrametric phylogenie(s) can have polytomie(s) (i.e. non binary tree).}
}
\value{
Returns the value of the stem age in the same unit as the branch length of the phylogenie(s) (typically in million of years).
}
\description{
Calculate the stem age of a rooted ultrametric phylogeny or set of phylogenies. This function can used only on stem crown phylogenie(s) and on binary or non binary phylogenie(s).
}
\details{
This function will calculate the stem age of a phylogeny or a set of phylogenies.
}
\seealso{
\code{\link{crown_age}}
}
\author{
Sophia Lambert
}
