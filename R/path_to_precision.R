#' Precision matrix from ordered path model
#'
#' Takes a path model and generates the corresponding precision matrix, assuming
#'    a standardised variance-covariance matrix. The inverse function to
#'    \code{\link{precision_to_path}}
#' @param B input \eqn{p \times p} weights matrix
#'
#' @return a \eqn{p \times p} precision matrix
#' @seealso \code{\link{precision_to_path}}, \code{\link{SEset_to_precision}}
#' @export
#' @importFrom Rdpack reprompt
#' @references
#'     \insertRef{ryan2019}{SEset}
#'
#'     \insertRef{shojaie2010penalized}{SEset}
#'
#'     \insertRef{bollen89sem}{SEset}
#' @examples

path_to_precision <- function(B){
  # first calculate lambda
  lambda <- solve(diag(nrow(B)) - B)
  # if(nrow(lambda)!=ncol(lambda)) return("Error: Lambda must be a square matrix")
  n <- nrow(lambda)
  psi <- diag(c(1,rep(0,n-1))) # Empty Matrix with one as first diagonal
  # Calculates psi such that a standardised sigma is outputted
  for(i in 2:n){
    psi[i,i] <- 1 - (t(lambda[i,]^2) %*% diag(psi))
  }
  # calculate sigma
  sigma <- lambda %*% psi %*% t(lambda)
  # calculate omega
  omega <- chol2inv(chol(sigma))
  dimnames(omega) <- dimnames(B)
  omega
}
