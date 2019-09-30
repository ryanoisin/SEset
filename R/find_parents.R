#' Return parent indices from a (weighted) DAG for a given child
#'
#' @param mat An \eqn{p \times p} weights or adjacency matrix
#' @param child Index giving the position of the child node
#' @return a vector containing index numbers defining the parent nodes
#' @seealso \code{\link{r2_distribution}}
#' @export
#' @references
#' \insertRef{ryan2019}{SEset}


find_parents <- function(mat, child){
    which(mat[child,] != 0)
  }
