#' Calculate scale for prior slab scales.
#'
#' The standard deviations for each column in the design matrix are to create
#' scaling values for the slab prior scale.
#'
#' @param x A design matrix.
#' @param min.x.sd The minimum acceptable standard deviation of each column.
#' @return A vector of scale values.
#' @importFrom stats sd
#' @note This function is taken from the R package \code{BhGLM}.
#' @return A vector of scaled values.
autoscale <- function (x, min.x.sd = 1e-04)
{
  scale <- apply(x, 2, sd)
  scale <- ifelse(scale <= min.x.sd, 1, scale)
  two <- which(apply(x, 2, function(u) {
    length(unique(u)) == 2
  }))
  scale[two] <- apply(x[, two, drop = FALSE], 2, function(u) {
    max(u) - min(u)
  })
  scale
}
