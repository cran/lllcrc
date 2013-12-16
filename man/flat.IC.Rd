\name{flat.IC}
\alias{flat.IC}
\title{Select an LLM}
\usage{
flat.IC(pop, rasch = FALSE, ic = "AICc", adjust = FALSE,
  averaging = FALSE)
}
\arguments{
  \item{pop}{A data.frame containing CRC data in standard
  form.}

  \item{rasch}{Logical: Should the Rasch model (most basic
  version, Darroch et. al. 1993) be considered, in addition
  to standard models?}

  \item{ic}{Character string specifying the information
  criterion to use for model selection.  Currently AIC,
  AICc, BIC, and BICpi are implemented.}

  \item{adjust}{Logical: Should we adjust the cells as in
  Evans and Bonett (1995)?}

  \item{averaging}{Logical: Should we use model averaging
  based on the information criterion scores?}
}
\value{
\item{pred}{The point estimate of the population size}
\item{form}{The log-linear terms in the chosen model}
}
\description{
Without using covariates (i.e., with capture probabilities
assumed flat over the covariate space), select the best
log-linear model for the marginal contingency table of
capture pattern counts.
}
\author{
Zach Kurtz
}
\references{
Fienberg SE (1972). "The Multiple Recapture Census for
Closed Populations and Incomplete $2^k$ Contingency
Tables." \emph{Biometrika}, \bold{59}(3), pp. 591.
}

