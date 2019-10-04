#' Path model from ordered precision matrix
#'
#' Takes a precision matrix and generates a lower-triangular weights matrix.
#'
#' @param omega input matrix, with rows and columns in desired topological ordering
#'     Must be an invertible square matrix
#' @param digits desired rounding of the output matrix
#' @param input_type specifies what type of matrix `omega` is.
#'     default is "precision", other options include a matrix of partial correlations
#'    ("parcor") or a model implied covariance or correlation matrix "MIcov"
#' @param quietly logical, show warning messages or not
#' @return lower triangular matrix containing regression weights of the path model.
#'   Element ij represents the effect of \eqn{X_j} on \eqn{X_i}
#' @seealso \code{\link{precision_to_SEset}}
#' @import stats
#' @export
#' @importFrom Rdpack reprompt
#' @references
#'     \insertRef{ryan2019}{SEset}
#'
#'     \insertRef{shojaie2010penalized}{SEset}
#'
#'     \insertRef{bollen89sem}{SEset}
#' @examples
#' data(riskcor)
#' omega <- (qgraph::EBICglasso(riskcor, n = 69, returnAllResults = TRUE))$optwi
#' # qgraph method estimates a non-symmetric omega matrix, but uses forceSymmetric to create
#' # a symmetric matrix (see qgraph:::EBICglassoCore line 65)
#' omega <- as.matrix(Matrix::forceSymmetric(omega)) # returns the precision matrix
#'
#' B <- precision_to_path(omega, digits=2)
#'
#' # Path model can be plotted as a weighted DAG
#' pos <- matrix(c(2,0,-2,-1,-2,1,0,2,0.5,0,0,-2),6,2,byrow=TRUE)
#'
#' # qgraph reads matrix elements as "from row to column"
#' # regression weights matrices are read "from column to row"
#' # path model weights matrix must be transposed for qgraph
#' qgraph::qgraph(t(B), labels=rownames(riskcor), layout=pos,
#' repulsion=.8, vsize=c(10,15), theme="colorblind", fade=FALSE)
#'


precision_to_path <- function(omega, input_type = "precision", digits=20, quietly = TRUE){

  # Save variable names
  label <- rownames(omega)

  # check that matrix is symmetric
  if(!Matrix::isSymmetric(omega)){
    omega <- as.matrix(Matrix::forceSymmetric(omega))
    if(!quietly){
    warning("Input matrix is not symmetric - Matrix::forceSymmetric() used to correct") }
    dimnames(omega) <- list(label,label)
  }

  # Take the inverse to produce a "model-implied" variance-covariance matrix
  if(input_type == "precision"){
  sigma <- chol2inv(chol(omega))
  diag(sigma) <- rep(1,dim(sigma)[1]) # set diagonal to exactly 1
  }
  # Else, transform the matrix of partial correlations
  if(input_type == "parcor"){
    if(!quietly){
    warning("If possible supply precision matrix as input")
    }
    omega <- -omega
    diag(omega) <- -diag(omega)
    sigma <- stats::cov2cor(solve(omega))
  }

  # Else, ensure that the supplied model implied covariance matrix is standardized
  if(input_type == "MIcov"){
     sigma <- cov2cor(omega)
  }

  # Solve for the LDL^T decomposition using the cholesky factor (GG^T)
  G <- chol(sigma)
  S <- diag(1/diag(G))
  U <- S %*% G  # This is the lambda ("influence") matrix

  # Solve for the regression weights matrix, rounded
  # backsolve to an adjacency matrix is quicker than solve with upper tri
  U <- round(diag(nrow(U)) - backsolve(U,x = diag(nrow(U))),digits = digits)
  if(!is.null(label)) {
    dimnames(U) <- list(label,label)} # Re-label
  # Matrix is in upper triangular form, transpose to lower
  t(U)

}
