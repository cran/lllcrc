% Generated by roxygen2 (4.0.2): do not edit by hand
\name{local.ml}
\alias{local.ml}
\title{Maximum likelihood estimation for fixed LLLMs}
\usage{
local.ml(pdat, ml.meth = NULL, predictors = NULL, k)
}
\arguments{
\item{pdat}{A design matrix with cell counts (possibly fractional)
included in the column named "c".}

\item{ml.meth}{A model specification such as "indep" or "indep.mono".}

\item{predictors}{A character vector of predictors of the form "c1", "c2"
for main effects, or "c12" for an interaction.}

\item{k}{Number of lists.}
}
\description{
A workhorse function for \code{apply.local.ml}.  For a fixed log-linear model,
get best fit using maximum local likelihood estimation.
}
\details{
Specify exactly one of the two arguments \code{predictors}, \code{ml.meth}
}
\author{
Zach Kurtz
}
\references{
Fienberg SE (1972). "The Multiple Recapture Census for Closed
Populations and Incomplete $2^k$ Contingency Tables." \emph{Biometrika},
\bold{59}(3), pp. 591.

Cormack RM (1989). "Log-linear models for capture-recapture."
\emph{Biometrics}, pp. 395-413.
}

