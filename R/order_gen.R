#' Generate all topological orderings
#'
#' Takes a matrix and generates a matrix containing all orderings of the rows and columns
#'
#' @param omega input p-dimensional square matrix
#' @return a \eqn{p \times p!} matrix of dimension orderings. Each column
#'     represents an ordering of dimension names as character strings.
#' @seealso \code{\link{reorder2}}, \code{\link{precision_to_SEset}}
#' @import combinat
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{combinat}{SEset}
#' @examples
#' data(riskcor)
#' orderings <- order_gen(riskcor)
#'
#' # Each column of orderings defines an ordering of variables
#' print(orderings[,1])
#' # in the second element, the fifth and sixth variable are switched
#' print(orderings[,2])


order_gen <- function(omega){
  if (dim(omega)[1] != dim(omega)[2]) {
    stop( "omega must be a square matrix" )
  }
  # if no dimension names are given, they are assigned here
  if (is.null(rownames(omega))) {
    dimnames(omega) <- list(paste0("X", seq( 1:nrow(omega))),
                            paste0("X",seq(1:nrow(omega))))
  }
  matrix(
    unlist( combinat::permn(rownames(omega)) )
    , nrow(omega), factorial(nrow(omega)))
}
