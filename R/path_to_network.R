#' Precision matrix from ordered path model
#'
#' Takes a path model and generates the corresponding (standardized) precision matrix or
#' covariance matrix. The inverse of \code{\link{network_to_path}}.
#' @param B input \eqn{p \times p} weights matrix
#' @param psi variance-covariance matrix for the residuals. If NULL (the default) will impose
#' the constraint that the variables have variance 1 and the residuals are uncorrelated
#' @param output Function returns the precision ("precision") or covariance ("covariance") matrix
#'
#' @return a \eqn{p \times p} precision or covariance matrix
#' @seealso \code{\link{network_to_path}}, \code{\link{SEset_to_network}}
#' @export
#' @importFrom Rdpack reprompt
#' @references
#'     \insertRef{ryan2019}{SEset}
#'
#'     \insertRef{shojaie2010penalized}{SEset}
#'
#'     \insertRef{bollen89sem}{SEset}

path_to_network <- function(B, psi = NULL, output = "precision"){
  # first calculate lambda
  lambda <- solve(diag(nrow(B)) - B)
  # if(nrow(lambda)!=ncol(lambda)) return("Error: Lambda must be a square matrix")
  n <- nrow(lambda)
  # if psi is not specified, impose constraint
  if(is.null(psi)){
    psi <- diag(c(1,rep(0,n-1))) # Empty Matrix with one as first diagonal
  # Calculates psi such that a standardised sigma is outputted
    for(i in 2:n){
    psi[i,i] <- 1 - (t(lambda[i,]^2) %*% diag(psi))
   }
  }
  # calculate sigma

  sigma <- lambda %*% psi %*% t(lambda)
  if(output == "covariance"){
  return(sigma)
  }else if(output == "precision"){
  # calculate omega
  omega <- chol2inv(chol(sigma))
  dimnames(omega) <- dimnames(B)
  return(omega)
  }
}
