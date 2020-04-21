#' Update conditional expectation of prior probabilities of inclusion at each iteration.
#'
#' @param group.vars A vector containing grouping information.
#' @param p The current estimate of conditional prior probabilities of inclusion.
#' @param p.bound A vector defining the lower and upper boundaries for the probabilities of inclusion
#' in the model, respectively. Defaults to \code{c(0.01, 0.99)}.
#' @return A list whose first element is a vector of updated scale parameters for each
#' parameter and whose second element is a vector of updated conditional expectations
#' of prior probabilities of model inclusion.
#' @note This function is taken unchanged from the R package \code{BhGLM}, except that bounds
#' for theta are now able to be specfied by the user rather than assumed to be 0.01 and 0.99.
#' @return A single numeric value between 0 and 1 that estimates the overall inclusion probability.
update_ptheta_group2 <- function(group.vars, p, p.bound = c(0.01, 0.99))
{
  theta <- p
  for (j in 1:length(group.vars)) {
    vars <- group.vars[[j]]
    theta[vars] <- mean(p[vars])
  }
  theta <- ifelse(theta < p.bound[1], p.bound[1], theta)
  theta <- ifelse(theta > p.bound[2], p.bound[2], theta)
  theta
}
