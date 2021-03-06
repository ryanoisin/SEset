% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_to_precision.R
\name{path_to_precision}
\alias{path_to_precision}
\title{Precision matrix from ordered path model}
\usage{
path_to_precision(B)
}
\arguments{
\item{B}{input \eqn{p \times p} weights matrix}
}
\value{
a \eqn{p \times p} precision matrix
}
\description{
Takes a path model and generates the corresponding precision matrix, assuming
   a standardised variance-covariance matrix. The inverse function to
   \code{\link{precision_to_path}}
}
\references{
\insertRef{ryan2019}{SEset}

    \insertRef{shojaie2010penalized}{SEset}

    \insertRef{bollen89sem}{SEset}
}
\seealso{
\code{\link{precision_to_path}}, \code{\link{SEset_to_precision}}
}
