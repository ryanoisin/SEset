#' Compute Controllability Distribution in the SE-set
#'
#' A function used to analyse the SEset results. For each member of the SE-set,
#' calculate the proportion of explained variance in each child node, when
#' predicted by all of its parent nodes
#'
#' @param SEmatrix An \eqn{n \times p} matrix containing the SEset, where each row
#'   represents a \eqn{p \times p} weights matrix stacked column-wise
#' @param cormat A \eqn{p \times p} matrix containing the marginal covariances or correlations
#' @param names optional character vector containing dimension names
#' @param indices option vector of matrix indices, indicating which variables to compute the R^2 distribution for
#'
#' @return Returns an \eqn{n \times p} matrix of \eqn{R^2} values.
#'  For each member of the SE-set, this represents the variance explained in node \eqn{X_i} by it's parents
#'  in that weighted DAG.
#' @seealso \code{\link{network_to_SEset}, \link{find_parents}}
#' @export
#' @references
#'     \insertRef{ryan2019}{SEset}
#'     \insertRef{haslbeck2018well}{SEset}
#' @examples
#' # first estimate the precision matrix
#' data(riskcor)
#' omega <- (qgraph::EBICglasso(riskcor, n = 69, returnAllResults = TRUE))$optwi
#' # qgraph method estimates a non-symmetric omega matrix, but uses forceSymmetric to create
#' # a symmetric matrix (see qgraph:::EBICglassoCore line 65)
#' omega <- as.matrix(Matrix::forceSymmetric(omega)) # returns the precision matrix
#'
#' SEmatrix <- network_to_SEset(omega, digits=3)
#'
#' r2set  <- r2_distribution(SEmatrix, cormat = riskcor, names = NULL, indices = c(1,3,4,5,6))
#' # Plot results
#' apply(r2set,2,hist)
#' # For ggplot format, execute
#' # r2set <- tidyr::gather(r2set)


r2_distribution <- function(SEmatrix, cormat, names = NULL, indices = NULL){
  # From a given correlation matrix and set of parents, find the R^2
  r2pa <- function(cormat = cormat, parents, child){
    # take the correlation of the parents and child
    Q <- t(cormat[child,parents])
    # take the inverse of parent submatrix
    subinv <- chol2inv(chol(cormat[parents,parents]))
    Q %*% subinv %*% t(Q)
  }

  # From an SE-set, find the R^2 over all DAGs, for a given child node
  r2find <- function(SEmatrix = SEmatrix, cormat = cormat, ind){
    p <- dim(cormat)[1]
    apply(SEmatrix, 1, function(row) {
      mat <- matrix(row,p,p)
      if (length(find_parents(mat, ind)) == 0) {
        return(0)
      } else {
        r2pa(cormat, parents = find_parents(mat, ind), child = ind)
      }
    })
  }

  # if indices is NULL, compute for all nodes
  if(is.null(indices)){ indices <- seq(1,sqrt(ncol(SEmatrix))) }
  # if names is NULL, create variable names
  if(is.null(names)){
    names <- colnames(cormat)[indices]
    if(is.null(names)){
    names <- paste0("X",indices)
    }}

  # Compute R^2 distribution
  out <- sapply(indices, function(x){ r2find(SEmatrix, cormat, ind = x) })
  dimnames(out) <- list(NULL, names)
  return(out)
}

