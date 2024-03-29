% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_to_path.R
\name{cov_to_path}
\alias{cov_to_path}
\title{Path model from covariance matrix with ordering}
\usage{
cov_to_path(sigma, ordering = NULL, digits = 2)
}
\arguments{
\item{sigma}{input matrix, with rows and columns in desired topological ordering
Must be an invertible square matrix}

\item{ordering}{character vector containing the dimension names of the input matrix
in the desired ordering}

\item{digits}{the number of digits used to round the output}
}
\value{
lower triangular matrix containing regression weights of the path model.
  Element ij represents the effect of \eqn{X_j} on \eqn{X_i}
}
\description{
Helper function. Takes a covariance matrix and ordering and generates a lower-triangular weights matrix.
}
\seealso{
\code{\link{network_to_path}}
}
