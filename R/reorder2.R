#' Re-order rows and columns
#'
#' Takes a matrix and re-orders the rows and columns to some target ordering
#'
#' @param matrix input matrix to be re-arranged. Must have rows and columns named
#' @param names character vector containing the dimension names of the input matrix
#'     in the desired ordering

#' @return input matrix with rows and columns sorted according to names
#' @seealso \code{\link{order_gen}}, \code{\link{precision_to_SEset}}
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{memisc}{SEset}
#' @examples
#' data(riskcor)
#'
#' # first define an ordered vector of names
#' names <- rownames(riskcor)
#' names <- c(names[1],names[2],names[3],names[4],names[6],names[5])
#'
#' reorder2(riskcor,names)
#'
#' # The fifth and sixth row and column have been switched
#' print(riskcor)

# function to reorder the matrix according to "names"
reorder2 <- function(matrix,names) {
  if (is.null(dimnames(matrix))) {
    stop("Error: matrix must have dimension names")
  }
  if (!all(rownames(matrix) %in% names)) {
    stop("Error: dimnames(matrix) does not match names")
  }
  # First re-order rows
  matrix_r <- memisc:::reorder.array(matrix,dim = 1,names = names)
  # Second re-order columns
  out <- memisc:::reorder.array(matrix_r, dim = 2,names = names)
  dimnames(out) <- list(names, names)
  out
}
