#' Calculate residual-covariance matrix based on a path model and covariance matrix
#'
#' Takes an ordered path model and corresponding variance-covariance matrix and computes the appropriate
#' residual covariance matrix (psi)
#' @param B input \eqn{p \times p} linear SEM weights matrix
#' @param sigma variance-covariance matrix of the variables
#' @param digits how many digits to round the result to
#'
#' @return a \eqn{p \times p} residual variance-covariance matrix
#' @export
#' @importFrom Rdpack reprompt

get_psi <- function(B, sigma, digits = 3){
  lambda <- solve(diag(nrow(B)) - B)
  psi  <- round(solve(lambda)%*%sigma%*%solve(t(lambda)),digits = digits)
  psi
}
