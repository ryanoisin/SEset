#' Path model from covariance matrix with ordering
#'
#' Helper function. Takes a covariance matrix and ordering and generates a lower-triangular weights matrix.
#'
#' @param sigma input matrix, with rows and columns in desired topological ordering
#'     Must be an invertible square matrix
#' @param ordering character vector containing the dimension names of the input matrix
#'     in the desired ordering
#' @param digits the number of digits used to round the output
#' @return lower triangular matrix containing regression weights of the path model.
#'   Element ij represents the effect of \eqn{X_j} on \eqn{X_i}
#' @seealso \code{\link{network_to_path}}
#' @export
#' @importFrom Rdpack reprompt

cov_to_path <- function(sigma, ordering = NULL, digits = 2){
  if(!is.null(ordering)){
    sigma <- SEset::reorder_mat(sigma,names = ordering)
  }
  # Solve for the LDL^T decomposition using the cholesky factor (GG^T)
  G <- chol(sigma)
  S <- diag(1/diag(G))
  U <- S %*% G  # This is the lambda ("influence") matrix

  # Solve for the regression weights matrix, rounded
  # backsolve to an adjacency matrix is quicker than solve with upper tri
  U <- round(diag(nrow(U)) - backsolve(U,x = diag(nrow(U))),digits = digits)
  out <- t(U)
  dimnames(out) = dimnames(sigma)
  out
}
