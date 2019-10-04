#' Precision matrices from the SEset
#'
#' Takes the SE-set and calculates for each weights matrix the corresponding
#'     precision matrix. Used to check the results of \code{\link{precision_to_SEset}}
#'     to assess deviations from statistical equivalence induced due to rounding,
#'     thresholding, and numerical approximations.
#' @param SEmatrix a \eqn{n \times p} matrix containing the SE-set. The output of
#'     \code{\link{precision_to_SEset}}
#' @param order.ref an optional character vector with variable names, the reference ordering
#'     of the precision matrix.
#' @param order.mat a \eqn{n \times p} matrix of character strings,
#'     defining the ordering of the matrix corresponding to each row of SEmatrix.
#'     If NULL it is assumed that all orderings are included and they are generated using
#'     \code{\link{order_gen}}
#' @param output Output as \code{"raw"} or \code{"summary"}. See value below
#' @param omega Comparision precision matrix, e.g. original input precision matrix to
#'     \code{\link{precision_to_SEset}}. Only necessary if \code{output = "summary"}
#' @return If \code{output = "raw"}, a \eqn{n \times p} matrix of precision matrices
#'     stacked column-wise in \eqn{n} rows.
#'     If \code{output = "summary"} returns a list containing the bias, MSE and
#'     RMSE for each re-calculated precision matrix, relative to comparison \code{omega}
#'     matrix supplied.
#' @seealso \code{\link{precision_to_path}}, \code{\link{path_to_precision}}
#' @export
#' @importFrom Rdpack reprompt
#' @references
#'     \insertRef{ryan2019}{SEset}
#'
#'     \insertRef{shojaie2010penalized}{SEset}
#'
#'     \insertRef{bollen89sem}{SEset}

SEset_to_precision <-
  function(SEmatrix, order.ref=NULL, order.mat = NULL, output="raw", omega=NULL){
  prec_bias <- function(precmat,omega) {

    comp <- omega[lower.tri(omega,diag = TRUE)]

    resid <- t(apply(precmat,1,function(row){row - comp}))
    sqresid <- t(apply(precmat,1,function(row) {(row - comp)^2} ))

    # Small function to make output look nicer
    matmake <- function(omega,residual){
      mat <- matrix(0,nrow(omega),nrow(omega))
      mat[lower.tri(mat,diag = TRUE)] <- colMeans(residual)
      mat[upper.tri(mat)] <- t(mat)[upper.tri(t(mat))]
      dimnames(mat) <- dimnames(omega)
      mat
    }

    bias <- matmake(omega,resid)
    MSE <- matmake(omega,sqresid)
    RMSE <- sqrt(MSE)

    results <- list(bias,MSE,RMSE)
    names(results) <- c("Bias","MSE","RMSE")
    results
  }


   nv <- sqrt(ncol(SEmatrix))
  # If no variable names are supplied, create some
  if (is.null(order.ref)) {
    order.ref <- paste0("X",seq(1:nv))
  }
  if (is.null(order.mat)) {
    order.mat <- t(matrix(unlist(combinat::permn(order.ref)), nv, factorial(nv)))
  }
  if (nrow(order.mat) != nrow(SEmatrix)) {
    return("Error: There must be one ordering per matrix (row) in SEmatrix")
  }

# Only the lower triangle & diagonal is necessary to save, symmetric matrices
  precmat <- matrix(0,nrow(order.mat),(nv^2 + nv)/2)
  for (i in 1:nrow(order.mat)) {
    Ab <- t(matrix(SEmatrix[i,],nv,nv))
    dimnames(Ab) <- list(order.ref,order.ref)
    Ab <- reorder_mat(Ab,names = order.mat[i,]) # matrix now topologically ordered
    omegaest <- reorder_mat(path_to_precision(Ab),order.ref)
    precmat[i,] <- omegaest[lower.tri(omegaest,diag = TRUE)]
  }
  if (output == "raw") {
    return(precmat)
  }else if (output == "performance") {
    if (is.null(omega)) {
      return("Error: for performance output, please supply the original
             precision matrix for comparison using the omega argument")
    } else {
    return(prec_bias(precmat,omega)) }
  }
}

