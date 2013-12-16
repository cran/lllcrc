\name{make.hierarchical.term.sets}
\alias{make.hierarchical.term.sets}
\title{Generate a universe of hierarchical model specifications.}
\usage{
make.hierarchical.term.sets(k = 3, rasch = FALSE)
}
\arguments{
  \item{k}{The number of lists}

  \item{rasch}{If TRUE, include the Rasch model (simplest
  Rasch model in Darroch et. al. 1993).}
}
\value{
A list of character vectors
}
\description{
Creates a list of all hierarchical sets of log-linear
modelling terms that include all main effects and exclude
the highest order interaction. Currently implemented only
for k=3 lists.
}
\details{
Needs generalization to k>3 very much!
}
\author{
Zach Kurtz
}

