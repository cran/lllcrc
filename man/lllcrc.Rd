\name{lllcrc}
\alias{lllcrc}
\title{Local log-linear models (LLLMs) for capture-recapture (CRC)}
\usage{
lllcrc(dat, kfrac, ic = "BICpi", bw = NULL, averaging = FALSE,
  round.vars = NULL, rounding.scale = 0.01, boot.control = NULL,
  cell.adj = TRUE)
}
\arguments{
  \item{dat}{Capture-recapture data, as output of
  \code{format.data}}

  \item{kfrac}{The approximate fraction of the data that is
  included in the support of the kernel for the local
  averages.}

  \item{ic}{The information criterion for selection of
  local log-linear models.  The default, BICpi, is
  introduced in and mentioned in Hook and Regal 1997.}

  \item{bw}{A matrix a single column, with rownames that
  match the covariate names in \code{dat}.  The values in
  the column are scalars that are used in constructing
  distances between covariate vectors.  Raw differences are
  divided by the corresponding scalars before being squared
  in the context of a Euclidean metric.}

  \item{averaging}{Logical: Should model averaging be done
  for each local model?}

  \item{round.vars}{See \code{micro.post.stratify}, which
  is called within \code{lllcrc}.}

  \item{rounding.scale}{See \code{micro.post.stratify},
  which is called within \code{lllcrc}.}

  \item{boot.control}{A list of control parameters for
  bootstrapping the sampling distribution of the
  estimator(s).  By default, there is no bootstrapping.}

  \item{cell.adj}{Logical: Should we adjust the cells as in
  Evans and Bonett (1995)?}
}
\value{
\item{est}{A point estimate of the population size}
\item{llform}{The set of log-linear terms} \item{dat}{The
output of function \code{micro.post.stratify}, with
estimated local rates of missingness appended as an extra
column labeled \code{pi0}.  In addition, \code{mct}
(multinomial cell count) gives the number of observed units
with that distinct covariate vector, and \code{cpi0}
(cumulative number missing) gives the the product of
\code{pi0} with \code{mct}, such that summing over this
vectorized product is exactly the Horvitz-Thompson style
sum in capture recapture. } \item{ess}{The local effective
sample sizes that are based on the local averaging weights
and used as eta_i in local model selection} \item{hpi}{The
matrix of local averages} \item{...}{The output is of class
\code{lllcrc} and has attributes \code{cont.x} and
\code{conteg.x}, which relate the continuous and
categorical variables in the model }
}
\description{
Fits local log-linear models.  Each distinct covariate
vector gets its own model.  To reduce the number of models,
some rounding of continuous covariates is done first.
}
\details{
The key implementation of the thesis of Kurtz 2013,
Carnegie Mellon University
}
\author{
Zach Kurtz
}
\references{
Kurtz ZT (2013). "Smooth Post-Stratification for Multiple
Capture-Recapture." \emph{arXiv preprint arXiv:1302.0890}.

Anderson DR and Burnham KP (1999). "Understanding
information criteria for selection among capture-recapture
or ring recovery models." \emph{Bird Study}, \bold{46}(S1),
pp. S14-S21.

Fienberg SE (1972). "The Multiple Recapture Census for
Closed Populations and Incomplete $2^k$ Contingency
Tables." \emph{Biometrika}, \bold{59}(3), pp. 591.
}

