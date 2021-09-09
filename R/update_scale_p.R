#' Update conditional expectation of prior probabilities of inclusion at each iteration.
#'
#' @param b0 A vector of current parameter estimates.
#' @param ss A vector of spike and slab prior scales, respectively.
#' @param theta The current estimate of prior probabilities of inclusion.
#' @param alpha A scalar value between 0 and 1 determining the compromise between the Ridge and Lasso models. When
#' \code{alpha = 1} reduces to the Lasso, and when \code{alpha = 0} reduces to Ridge.
#' @return A list whose first element is a vector of updated scale parameters for each
#' parameter and whose second element is a vector of updated conditional expectations
#' of prior probabilities of model inclusion.
#' @note This function is a modified version of \code{update_scale_p()} from the R package \code{BhGLM}.
update_scale_p <- function (b0, ss, theta, alpha)
{
  if (alpha == 1) {
    den0 <- (2 * ss[1])^(-1) * exp(-abs(b0)/ss[1])
    den1 <- (2 * ss[2])^(-1) * exp(-abs(b0)/ss[2])
  } else {
    den0 <- (1 - alpha) * exp(-log(sqrt(2 * pi * ss[1])) - b0^2 / ss[1]) +
            alpha * exp(-log(2 * ss[1]) - abs(b0) / ss[1])
    den1 <- (1 - alpha) * exp(-log(sqrt(2 * pi * ss[2])) - b0^2 / ss[2]) +
      alpha * exp(-log(2 * ss[2]) - abs(b0) / ss[2])
  }
  p <- theta * den1/(theta * den1 + (1 - theta) * den0 + 1e-10)
  scale <- 1/((1 - p)/ss[1] + p/ss[2] + 1e-10)
  list(scale = scale, p = p)
}
