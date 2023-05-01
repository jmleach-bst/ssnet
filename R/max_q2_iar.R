#' Update prior probabilities of model inclusion at each iteration.
#'
#' Uses \code{stan} to optimize the term in the EM algorithm that depends on
#' the prior probabilities of inclusion in the model, \eqn{p_j}.
#'
#' @param iar.data A list of output from \code{\link{mungeCARdata4stan}} that
#' contains the necessary inputs for the IAR prior.
#' @param p A vector containing the prior probabilities of inclusion for each
#' parameter at the current iteration.
#' @param opt.algorithm One of the optimization algorithms available from
#' \code{optimizing}.
#' @param tau.prior One of \code{c("none", "manual", "cauchy")}. This argument
#' determines the precision parameter in the Conditional Autoregressive model
#' for the (logit of) prior inclusion probabilities. When \code{"none"}, the
#' precision is set to 1; when "manual", the precision is manually entered by
#' the user (and still not random); when \code{"cauchy"}, the inverse precision
#' is assumed to follow a Cauchy distribution with mean 0 and scale 2.5. Note
#' that at this stage of development, only the \code{"none"} option has been
#' extensively tested, so the other options should be used with caution.
#' @param tau.manual When \code{tau.prior = "manual"}, use this argument to
#' specify a common precision parameter.
#' @param p.bound A vector defining the lower and upper boundaries for the
#' probabilities of inclusion in the model, respectively. Defaults to
#' \code{c(0.01, 0.99)}.
#' @param stan_manual A \code{stan_model} that is manually specified.
#' @importFrom Rdpack reprompt
#' @return A vector containing updated probabilities of inclusion.
#' @note This function borrows from the work of Mitzi Morris, who describes
#' how to fit an intrinsic autoregression in
#' \code{stan} \insertCite{Morris:2017,Morris:2019}{ssnet}.
#' @references
#'
#' \insertRef{Morris:2017}{ssnet}
#'
#' \insertRef{Morris:2019}{ssnet}
#'
#' \insertRef{Tang:2017}{ssnet}
max_q2_iar <- function(
    iar.data,
    p,
    opt.algorithm = "LBFGS",
    tau.prior = "none",
    tau.manual = NULL,
    p.bound = c(0.01, 0.99),
    stan_manual = NULL
){
  p[p < p.bound[1]] <- p.bound[1]
  p[p > p.bound[2]] <- p.bound[2]
  iar.data$p <- p
  if (is.null(stan_manual) == FALSE) {
    Q2 <- rstan::optimizing(
      object = stan_manual,
      data = iar.data,
      algorithm = opt.algorithm
    )
  } else {
      if (tau.prior == "cauchy") {
        Q2 <- rstan::optimizing(
          object = stanmodels$iar_incl_prob_current,
          data = iar.data,
          algorithm = opt.algorithm
        )
      }
      if (tau.prior == "manual") {
        iar.data$tau <- tau.manual
        Q2 <- rstan::optimizing(
          object = stanmodels$iar_incl_prob_manual_tau,
          data = iar.data,
          algorithm = opt.algorithm
        )
      }
      if (tau.prior == "none") {
        Q2 <- rstan::optimizing(
          object = stanmodels$iar_incl_prob_notau,
          data = iar.data,
          algorithm = opt.algorithm
        )
      }
    }
  theta <- Q2$par[grep("theta*", names(Q2$par))]
  theta[theta < p.bound[1]] <- p.bound[1]
  theta[theta > p.bound[2]] <- p.bound[2]
  return(theta)
}
