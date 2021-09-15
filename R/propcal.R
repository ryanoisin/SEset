#' Edge frequency in the SE-set
#'
#' A function used to analyse the SEset results. Calculates the proportion of
#'     path models in a given SEset in which a particular edge is present
#'
#' @param SEmatrix An \eqn{n \times p} matrix containing the SEset, where each row
#'   represents a \eqn{p \times p} weights matrix stacked column-wise
#' @param names optional character vector containing dimension names
#' @param rm_duplicate Should duplicate weights matrices be removed from the SEset.
#'     Defaults to TRUE.
#' @param directed If \code{FALSE}, the directionality of edges is ignored, and the output
#'    reflects in what proportion of the SEset an edge of any direction is present.
#'    If \code{TRUE}, the proportion is calculated seperately for edges of either direction.
#'    Defaults to TRUE
#' @return a \eqn{p \times p} matrix showing in what proportion particular edges are present.
#'    If directed=TRUE, elements ij denote the proportion of weights matrices containing a path
#'    from \eqn{X_j} to \eqn{X_i}. If directed=F, the output will be a symmetric matrix, with element ij
#'    denoting in what proprtion an edge of either direction connects \eqn{X_i} to \eqn{X_j}.
#' @seealso \code{\link{network_to_SEset}}
#' @export
#' @references
#' \insertRef{ryan2019}{SEset}


propcal <- function(SEmatrix, names=NULL, rm_duplicate = TRUE, directed = TRUE){
  if (is.null(dim(SEmatrix))) {
    stop("SEset must have more than one element")
  }

  if (rm_duplicate == TRUE) {
    SEmatrix <- unique(SEmatrix,MARGIN = 1)
  }
  SEmatrix[SEmatrix != 0] <- 1
  n <- sqrt(ncol(SEmatrix))
  # matrix is in Upper-triangular form (transpose of SEM beta matrix)
  prop <- matrix(apply(SEmatrix,2,sum)/dim(SEmatrix)[1],n,n)
  if (directed == F) {
    prop = prop + t(prop)
  }
  if (!is.null(names)) {
    dimnames(prop) <- list(names,names)
  }
  prop
  #lower triangular form - transpose for qgraph
}

