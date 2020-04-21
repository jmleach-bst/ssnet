#' Update conditional expectation of prior probabilities of inclusion at each iteration.
#'
#' @param b0 A vector of current parameter estimates.
#' @param ss A vector of spike and slab prior scales, respectively.
#' @param theta The current estimate of prior probabilities of inclusion.
#' @return A list whose first element is a vector of updated scale parameters for each
#' parameter and whose second element is a vector of updated conditional expectations
#' of prior probabilities of model inclusion.
#' @note This function is taken unchanged (except for the name) from the R package \code{BhGLM}.
update_scale_p <- function (b0, ss, theta)
{
  den0 <- (2 * ss[1])^(-1) * exp(-abs(b0)/ss[1])
  den1 <- (2 * ss[2])^(-1) * exp(-abs(b0)/ss[2])
  p <- theta * den1/(theta * den1 + (1 - theta) * den0 + 1e-10)
  scale <- 1/((1 - p)/ss[1] + p/ss[2] + 1e-10)
  list(scale = scale, p = p)
}
