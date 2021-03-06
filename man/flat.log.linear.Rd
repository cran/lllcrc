% Generated by roxygen2 (4.0.2): do not edit by hand
\name{flat.log.linear}
\alias{flat.log.linear}
\title{Fit an LLM}
\usage{
flat.log.linear(pop, model.terms, rasch = FALSE)
}
\arguments{
\item{pop}{The CRC data as a data frame.}

\item{model.terms}{The columns of the standard design matrix to include in
the model.  For example, "c1", "c2" for main effects, and "c12" for
interactions.}

\item{rasch}{Logical: Is this the Rasch model?}
}
\value{
A vector of log-linear coefficients.  The first coefficient is the
intercept, and the rest correspond (in order) with the \code{model.terms}
argument
}
\description{
Fit a log-linear model.  This is a wrapper function for our own variant of
the \code{glm} function, \code{\link{pirls}}.
}
\details{
Maximum likelihood estimation is used, conditioning on the observed
population as if it were the full population.
}
\author{
Zach Kurtz
}

