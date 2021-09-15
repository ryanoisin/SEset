#' Re-order rows and columns
#'
#' Takes a matrix and re-orders the rows and columns to some target ordering
#'
#' @param matrix input matrix to be re-arranged. Must have rows and columns named
#' @param names character vector containing the dimension names of the input matrix
#'     in the desired ordering

#' @return input matrix with rows and columns sorted according to names
#' @seealso \code{\link{order_gen}}, \code{\link{network_to_SEset}}
#' @export
#' @importFrom Rdpack reprompt
#' @examples
#' data(riskcor)
#'
#' # first define an ordered vector of names
#' row_names <- rownames(riskcor)
#' row_names_new <- row_names[c(1,2,3,4,6,5)]
#'
#' reorder_mat(riskcor,row_names_new)
#'
#' # The fifth and sixth row and column have been switched
#' print(riskcor)

# function to reorder the matrix according to "names"
reorder_mat <- function(matrix, names) {
  if (is.null(dimnames(matrix))) {
    stop("Error: matrix must have dimension names")
  }
  if (!all(rownames(matrix) %in% names)) {
    stop("Error: dimnames(matrix) does not match names")
  }
  matrix[names,names]
}
