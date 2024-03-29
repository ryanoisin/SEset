#' SE-set from precision matrix
#'
#' Takes a precision matrix and generates the SE-set, a set of
#'     statistically equivalent path models. Unless otherwise specified, the SEset
#'     will contain one weights matrix for every possible topological ordering of the
#'    input precision matrix
#' @param omega input \eqn{p \times p} precision matrix
#' @param orderings An optional matrix of \eqn{n} orderings from which to generate the
#'     SE-set. Must be in the form of a \eqn{p \times n} matrix with each column a vector of
#'     dimension names in the desired order. If unspecified, all \eqn{p!} possible
#'     orderings are used
#' @param digits desired rounding of the output weights matrices in the SE-set,
#'     in decimal places. Defaults to 20.
#' @param rm_duplicates Logical indicating whether only unique DAGs should be returned
#' @param input_type specifies what type of matrix `omega` is.
#'     default is "precision", other options include a matrix of partial correlations
#'    ("parcor") or a model implied covariance or correlation matrix "MIcov"
#' @return a \eqn{p! \times p} matrix containing the SE-set
#'     (or \eqn{n \times p}  matrix if a custom set of \eqn{n} orderings is specified).
#'     Each row represents a lower-triangular weights matrix, stacked column-wise.
#' @seealso \code{\link{network_to_path}}, \code{\link{reorder_mat}}, \code{\link{order_gen}}
#' @export
#' @importFrom Rdpack reprompt
#' @references
#'     \insertRef{ryan2019}{SEset}
#'
#'     \insertRef{shojaie2010penalized}{SEset}
#'
#'     \insertRef{bollen89sem}{SEset}
#' @examples
#' # first estimate the precision matrix
#' data(riskcor)
#' omega <- (qgraph::EBICglasso(riskcor, n = 69, returnAllResults = TRUE))$optwi
#' # qgraph method estimates a non-symmetric omega matrix, but uses forceSymmetric to create
#' # a symmetric matrix (see qgraph:::EBICglassoCore line 65)
#' omega <- as.matrix(Matrix::forceSymmetric(omega)) # returns the precision matrix
#'
#' SE <- network_to_SEset(omega, digits=3)
#'
#' # each row of SE defines a path-model weights matrix.
#' # We can extract element 20 by writing it to a matrix
#' example <- matrix(SE[20,],6,6)
#'
#' # Example path model can be plotted as a weighted DAG
#' pos <- matrix(c(2,0,-2,-1,-2,1,0,2,0.5,0,0,-2),6,2,byrow=TRUE)
#'
#' # qgraph reads matrix elements as "from row to column"
#' # regression weights matrices are read "from column to row"
#' # path model weights matrix must be transposed for qgraph
#'
#' qgraph::qgraph(t(example), labels=rownames(riskcor), layout=pos,
#' repulsion=.8, vsize=c(10,15), theme="colorblind", fade=FALSE)

network_to_SEset <- function(omega, orderings=NULL, digits=20, rm_duplicates = FALSE,
                               input_type = "precision"){

  # check that matrix is symmetric
  if(!Matrix::isSymmetric(omega)){
    label <- dimnames(omega)
    omega <- as.matrix(Matrix::forceSymmetric(omega))
    warning("Input matrix is not symmetric - Matrix::forceSymmetric() used to correct")
    dimnames(omega) <- label
  }

  # If variables are unnamed, name them here

  if (is.null(dimnames(omega)[[1]])) {
    dimnames(omega) <- list(paste0("X",seq(1:nrow(omega))),
                            paste0("X",seq(1:nrow(omega))))
  }

  if(input_type == "parcor"){
  }
  # If no orderings supplied, all possible orderings are taken
  if (is.null(orderings)) {
    orderings <- order_gen(omega)
  }
  # For each ordering, calculate the adjcacency matrix of the DAG
  out <- t(apply(orderings,2,function(per) {
    omega_r <- reorder_mat(omega,names = per)
    reorder_mat((network_to_path(omega_r, input_type = input_type, digits = digits)),
             names = rownames(omega))
  }) )
  if (rm_duplicates){
    out <- unique(out, MARGIN = 1)
  }
  return(out)

}

