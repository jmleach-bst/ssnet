#' Update conditional expectation of prior probabilities of inclusion at each iteration.
#'
#' Specifically adresses group selection where conditional expecation depends on vector of parameters, not
#' a single scalar parameter value, i.e., this was specifically developed for our multinomial extension of
#' spike-and-slab EN GLMs.
#'
#' @param b0 A list. The length of the list should be equal to the total number of predictors, and each element
#' should be a vector of parameter estimates corresponding to each
#' @param ss A vector of spike and slab prior scales, respectively.
#' @param theta The current estimate of prior probabilities of inclusion.
#' @param alpha A scalar value between 0 and 1 determining the compromise between the Ridge and Lasso models. When
#' \code{alpha = 1} reduces to the Lasso, and when \code{alpha = 0} reduces to Ridge.
#' @return A list whose first element is a vector of updated scale parameters for each
#' parameter and whose second element is a vector of updated conditional expectations
#' of prior probabilities of model inclusion.
#' @note This function is a modified version of \code{update_scale_p()} from the R package \code{BhGLM}.
update_scale_p_mn <- function(b0, ss, theta, alpha, print_out = FALSE) {

  # b0 is initially a list with number of elements equal to number of outcome categories, V,
  # and each element of length J. We need to make this list into a matrix and transpose it,
  # so that each row contains a V-length vector of parameters for a total of J rows.
  b0 <- matrix(
    unlist(b0),
    ncol = length(b0),
    byrow = FALSE
  )

  if (print_out == TRUE) {
    cat("Current parameter estimates \n")
    print(b0)
  }

  if (nrow(b0) != length(theta)) {
    stop("nrow(b0) != length(theta)")
  }

  # DE prior density
  de.pdf <- function(x, s) {
    (1 / (2 * s)) * exp(-abs(x) / s)
  }

  # EN prior density
  en.pdf <- function(x, s, a) {
    (1 - alpha) * exp(-log(sqrt(2 * pi * s)) - (x^2 / s)) +
      alpha * exp(-log(2 * s) - abs(x) / s)
  }

  ##################################################
  # code for p_j in MN paper (i.e., equation (12)) #
  ##################################################

  if (alpha == 1) {
    # p(b_j^v | gamma_j = 0, s0)
    p0.b0 <- apply(
      X = b0, MARGIN = 2, FUN = de.pdf, s = ss[1]
    )

    # p(b_j^v | gamma_j = 1, s1)
    p1.b0 <- apply(
      X = b0, MARGIN = 2, FUN = de.pdf, s = ss[2]
    )

  } else {
    # p(b_j^v | gamma_j = 0, s0)
    p0.b0 <- apply(
      X = b0, MARGIN = 2, FUN = en.pdf, s = ss[1]
    )

    # p(b_j^v | gamma_j = 1, s1)
    p1.b0 <- apply(
      X = b0, MARGIN = 2, FUN = en.pdf, s = ss[2]
    )
  }

  if (print_out == TRUE) {
    cat("p(b_j^v | gamma_j = 0, s0) \n")
    print(p0.b0)
    cat("p(b_j^v | gamma_j = 1, s1) \n")
    print(p1.b0)
  }

  # prod_{v=1}^{V} p(b_j^v | gamma_j = 0, s0))
  prod.p0.b0 <- apply(
      X = p0.b0, MARGIN = 1, FUN = matrixStats::product
    )
  # p(gamma_j = 0 | theta_j) * prod_{v=1}^{V} p(b_j^v | gamma_j = 0, s0)
  pt.p0.b0 <- (1 - theta) * prod.p0.b0

  # prod_{v=1}^{V} p(b_j^v | gamma_j = 1, s1))
  prod.p1.b0 <- apply(
      X = p1.b0, MARGIN = 1, FUN = matrixStats::product
    )
  # p(gamma_j = 1 | theta_j) * prod_{v=1}^{V} p(b_j^v | gamma_j = 1, s0)
  pt.p1.b0 <- theta * prod.p1.b0

  if (print_out == TRUE) {
    cat("prod_{v=1}^{V} p(b_j^v | gamma_j = 0, s0)) \n")
    print(prod.p0.b0)
    cat("prod_{v=1}^{V} p(b_j^v | gamma_j = 1, s1)) \n")
    print(prod.p1.b0)
  }

  # conditional probabilities p_j = p(\gamma_j = 1 | \bm{\beta}_j, \theta_j)
  p <- pt.p1.b0 / (pt.p1.b0 + pt.p0.b0 + 1e-10)
  scale <- 1/((1 - p)/ss[1] + p/ss[2] + 1e-10)
  return(list(scale = scale, p = p))
}
