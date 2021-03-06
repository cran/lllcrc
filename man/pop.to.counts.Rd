% Generated by roxygen2 (4.0.2): do not edit by hand
\name{pop.to.counts}
\alias{pop.to.counts}
\title{Put CRC data into LLM vector}
\usage{
pop.to.counts(y)
}
\arguments{
\item{y}{A vector of capture patterns}
}
\description{
Essentially, this does \code{table(y)}, where \code{y} is a vector of
capture patterns, with the exception that here every capture pattern that is
observable (in principle) is included, with a possibly zero count.
}
\author{
Zach Kurtz
}

