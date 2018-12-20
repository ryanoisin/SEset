#' Estimate the precision matrix
#'
#' A function to compute a sparse gaussian graphical model. This function is
#'     taken directly from the \link[qgraph]{qgraph} package function \link[qgraph]{EBICglasso}
#'     version 1.4.4 (Epskamp et al. 2012).
#'     The only alteration is to remove the wi2net function call so as to output
#'     the precision matrix directly, rather than the standardised partial correlations
#'
#' @param S A covariance or correlation matrix
#' @param n Sample size used in computing \code{S}
#' @param gamma EBIC tuning parameter. 0.5 is generally a good choice.
#'     Setting to zero will cause regular BIC to be used.
#' @param penalize.diagonal Should the diagonal be penalized?
#' @param nlambda Number of lambda values to test.
#' @param lambda.min.ratio Ratio of lowest lambda value compared to maximal lambda
#' @param returnAllResults  If \code{TRUE} this function does not return a network
#'      but the results of the entire glasso path.
#' @param checkPD If \code{TRUE}, the function will check if \code{S} is
#'     positive definite and return an error if not. It is not advised to use
#'      a non-positive definite matrix as input as (a) that can not be a covariance
#'      matrix and (b) glasso can hang if the input is not positive definite.
#' @param penalizeMatrix  Optional logical matrix to indicate which elements are penalized
#' @param countDiagonal Should diagonal be counted in EBIC computation? Defaults to \code{FALSE}.
#'     Set to \code{TRUE} to mimic qgraph < 1.3 behavior (not recommended!).
#' @param refit Logical, should the optimal graph be refitted without LASSO regularization?
#'     Defaults to \code{FALSE}.
#' @param ... functions to be passed to \link[glasso]{glasso}
#' @return a \eqn{p \times p} estimated precision matrix
#' @seealso
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{qgraph}{SEset}
#'
#' \insertRef{glasso}{SEset}
#'
#' \insertRef{friedman2008sparse}{SEset}
#'
#' \insertRef{Matrix}{SEset}
#'
#' \insertRef{ryan2019}{SEset}
#'
#' \insertRef{hoorelbeke2016interplay}{SEset}
#'
#' @examples
#' data(riskcor)
#' omega <- EBICglasso_wi(riskcor,n=69)
#'
#' # Results can be plotted as a GGM
#' parcor <- wi2net(omega)
#' pos <- matrix(c(2,0,-2,-1,-2,1,0,2,0.5,0,0,-2),6,2,byrow=T)
#' qgraph(parcor, labels=rownames(riskcor), layout=pos,
#' repulsion=.8, vsize=c(10,15), theme="colorblind")

EBICglasso_wi <- function (S, n, gamma = 0.5, penalize.diagonal = FALSE, nlambda = 100,
          lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE,
          penalizeMatrix, countDiagonal = FALSE, refit = FALSE, ...)
{
  if(!require("Matrix")) install.packages("Matrix")
  if(!require("glasso")) install.pacakges("glasso")
  if(!require("qgraph")) install.packages("qgraph")

  require(Matrix)
  # require(qgraph)
  require(glasso)
  require(qgraph)

  ## According to huge???
  logGaus <- function(S,K,n)
  {
    KS = K %*% S
    # SK = S %*% K
    tr = function(A) sum(diag(A))
    return(n/2 * (log(det(K)) - tr(KS))  )
    # return(n/2 * (log(det(K)) - tr(SK))  )
  }

  # Computes the EBIC:
  EBIC <- function(S,K,n,gamma = 0.5,E,countDiagonal=FALSE)
  {
    #   browser()
    L <- logGaus(S, K, n)
    if (missing(E)){
      E <- sum(K[lower.tri(K,diag=countDiagonal)] != 0)
      # E <- sum(abs(K[lower.tri(K,diag=countDiagonal)]) > sqrt(.Machine$double.eps))
    }
    p <- nrow(K)

    # return EBIC:
    -2 * L + E * log(n) + 4 * E * gamma * log(p)
  }

  if (checkPD) {
    if (any(eigen(S)$values < 0))
      stop("'S' is not positive definite")
  }
  S <- cov2cor(S)
  lambda.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  lambda.min = lambda.min.ratio * lambda.max
  lambda = exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  if (missing(penalizeMatrix)) {
    glas_path <- glassopath(S, lambda, trace = 0, penalize.diagonal = penalize.diagonal,
                            ...)
  }
  else {
    glas_path <- list(w = array(0, c(ncol(S), ncol(S), length(lambda))),
                      wi = array(0, c(ncol(S), ncol(S), length(lambda))),
                      rholist = lambda)
    for (i in 1:nlambda) {
      res <- glasso(S, penalizeMatrix * lambda[i], trace = 0,
                    penalize.diagonal = penalize.diagonal, ...)
      glas_path$w[, , i] <- res$w
      glas_path$wi[, , i] <- res$wi
    }
  }
  lik <- sapply(seq_along(lambda), function(i) {
    logGaus(S, glas_path$wi[, , i], n)
  })
  EBICs <- sapply(seq_along(lambda), function(i) {
    EBIC(S, glas_path$wi[, , i], n, gamma, countDiagonal = countDiagonal)
  })
  opt <- which.min(EBICs)
  if (opt == 1) {
    warning("Network with lowest lambda selected as best network. Try setting 'lambda.min.ratio' lower.")
  }
  net <- as.matrix(forceSymmetric(glas_path$wi[, , opt]))
  colnames(net) <- rownames(net) <- colnames(S)
  if (all(net == 0)) {
    message("An empty network was selected to be the best fitting network. Possibly set 'lambda.min.ratio' higher to search more sparse networks. You can also change the 'gamma' parameter to improve sensitivity (at the cost of specificity).")
  }
  if (refit) {
    message("Refitting network without LASSO regularization")
    glassoRes <- suppressWarnings(glasso::glasso(S, 0, zero = which(net ==
                                                                      0 & upper.tri(net), arr.ind = TRUE), trace = 0, penalize.diagonal = penalize.diagonal,
                                                 ...))
    net <- as.matrix(forceSymmetric(glassoRes$wi))
    colnames(net) <- rownames(net) <- colnames(S)
    optwi <- glassoRes$wi
  }
  else {
    optwi <- glas_path$wi[, , opt]
  }
  if (returnAllResults) {
    return(list(results = glas_path, ebic = EBICs, loglik = lik,
                optnet = net, lambda = lambda, optwi = optwi))
  }
  else return(net)
}



