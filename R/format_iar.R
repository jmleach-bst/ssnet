#' Format IAR data for ssnet
#'
#' Given an adjacency matrix, formats data for using IAR priors to model spatial structure. Otherwise, assumes
#' default structure where neighbors are defined as variables directly above or to the side of a specified
#' variable location.
#' @param x Design, or input, matrix, of dimension nobs x nvars; each row is an observation vector. It is recommended that
#' \code{x} have user-defined column names for ease of identifying variables.
#' @param adjmat A data.frame or matrix containing a "sparse" representation of the neighbor relationships. The first
#' column should contain a numerical index for a given location. Each index will be repeated in this column for
#' every neighbor it has. The indices for the location's neighbors are then specified in the second column. An optional
#' third column specifies weights. If no third column is specified, then equal weights are assumed.
#' @param im.res A 2-element vector where the first argument is the number of "rows" and the second argument
#' is the number of "columns" in each subject's "image". Default is \code{NULL}.
#' @return A list containing data formatted for using IAR priors.
#'
format_iar <- function(adjmat = NULL, im.res = NULL, x = NULL,
                       tau.prior = "none", tau.manual = NULL) {
  # construct adjacency matrix
  if (is.null(adjmat) == TRUE) {
    if (is.null(im.res) == TRUE | is.null(x) == TRUE) {
      stop("User must specify both im.res and x.")
    } else {
      if (ncol(x) != prod(im.res)) {
        stop("The product of im.res should equal the number of columns in x.")
      }
    }
    adjmat <- sim2Dpredictr::proximity_builder(im.res = im.res, type = "sparse")
  } else {
    adjmat <- as.data.frame(adjmat)
    colnames(adjmat)[1:2] <- c("location.index", "nb.index")
  }

  if (tau.prior == "manual") {
    if (is.null(tau.manual) == TRUE) {
      stop("User must specify tau.manual")
    }
    mcd <- mungeCARdata4stan(adjmat$nb.index,
                             table(adjmat$location.index))
    mcd$tau <- tau.manual
    return(mcd)
  } else {
    mcd <- return(mungeCARdata4stan(adjmat$nb.index,
                                    table(adjmat$location.index)))
  }

}
