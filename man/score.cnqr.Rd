% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.cnqr.R
\name{score.cnqr}
\alias{score.cnqr}
\title{Score CNQR Models}
\usage{
\method{score}{cnqr}(obj, dat = "tr", sc = NULL, QY = NULL, w = 1,
  cmb = TRUE, se = FALSE, na_omit = TRUE)
}
\arguments{
\item{obj}{Object of type \code{'cnqr'}. See \code{\link{cnqr}}.}

\item{dat}{Data matrix or data frame with raw data.
Or, a character vector specifying data that's already in \code{'obj'} --
for example, \code{"tr"} for the training data (always present), or
\code{"val"} for validation data. The outputted score is of the
pooled data in the \code{dat} vector.}

\item{sc}{Scoring rule to use for the regression, as in the output
of \code{\link{scorer}}.}

\item{QY}{(Quantile) function that accepts a vector of quantile levels, and
returns an equally long vector of corresponding quantiles (only a
stationary quantile function is accepted at the moment)
Alternatively, \code{NULL} if the quantile function is
already in \code{obj} (and you specify \code{dat} as a character).}

\item{w}{Vector of weights corresponding to the observations in
\code{yhat}. No need to normalize them. See details to see exactly how these
weights are used.}

\item{cmb}{Logical; should the scores be combined (via average)? \code{TRUE}
if so, \code{FALSE} to output a score matrix for each observation (rows)
and each quantile level (column).}

\item{se}{Logical; should an estimate of the standard error of the mean score
estimate be returned? \code{TRUE} if so. Is not considered if
\code{cmb} is \code{FALSE}.}

\item{na_omit}{Logical; should observations leading to an \code{NA} score
(for any \code{tau}) be removed? \code{TRUE} by default. Warning message
appears when observations are removed.}
}
\value{
If \code{cmb} is \code{FALSE}, returns a matrix of scores with rows
corresponding to observations, and columns
corresponding to quantile levels.

If \code{cmb} is \code{TRUE}, then either a single numeric of the
average score is returned (if \code{se} is \code{FALSE}), or a named vector
of the mean score estimate, followed by an estimate of its standard error
(using the sample standard deviation through the \code{sd} function, divided
by sqrt of number of observations).
}
\description{
Assesses \code{"cnqr"} fits using a Proper Scoring Rule.
}
\examples{
data(egdat)
dat <- list(egdat[1:750, ], egdat[751:1000, ])
basevine <- subset(egvine, 1:4)
sc <- scorer(space_taus(10))
fit <- cnqr(5:4, dat, basevine, sc=sc)

score(fit)
score(fit, dat=egdat, QY=identity)
head(score(fit, cmb=FALSE))
}
