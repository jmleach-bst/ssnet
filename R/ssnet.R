#' Bayesian Spike-and-Slab Elastic Net with Spatial Structure
#'
#' Fits generalized linear models with spike-and-slab priors whose (logit of) probability of inclusion in the model
#' is either assigned an intrinsic autoregression as a prior to incorporate spatial information or is unstructured.
#' The model is fit using an EM algorithm where the E-step is fit by \code{glmnet()} and the M-step is fit using the
#' \code{stan} function \code{optimizing} (when IAR prior is employed).
#'
#' @param verbose Logical. If \code{TRUE}, prints out the number of iterations and computational time.
#' @importFrom Rdpack reprompt
#' @importFrom stats sd
#' @importFrom utils install.packages
#' @return An object of class \code{c("elnet"   "glmnet"  "bmlasso" "GLM")}.
#' @note If \code{iar.data = NULL}, i.e. is left unspecified, then provided that \code{im.res} is specified, the function
#' \code{proximity_builder()} from the package \code{sim2Dpredictr} builds the appropriate list of data for
#' optimization with \code{stan}. Currently, \code{im.res} can only handle 2D data. Future versions may allow
#' images to be 3D. However, the function will work given any appropriately specified neighborhood matrix,
#' whatever the original dimension.
#' @examples
#' library(sim2Dpredictr)
#' ## build adjacency matrix
#' adjmat10x10 <- sim2Dpredictr::proximity_builder(im.res = c(10, 10), type = "sparse")
#'
#' ## stan model information
#' model_info10x10 <- mungeCARdata4stan(adjmat10x10$nb.index,
#'                                      table(adjmat10x10$location.index))
#' ## pre-specify stan model
#' sm <- stan_model(file =
#' "C:/Users/Justin/Documents/BST/Dissertation_in_Latex/stan models/iar_incl_prob_notau.stan"
#' )
#'
#' ## fit model
#' cn <- c()
#' for (i in 1:100) cn[i] <- paste0("x", i)
#' bmn_model <- ssnet(x = matrix(rnorm(10000), nrow = 100, ncol = 100,
#'                    dimnames = list(1:100, cn)), alpha = 0.5,
#'                    iar.prior = TRUE,
#'                    y = rnorm(100), iar.data = model_info10x10,
#'                    family = "gaussian", stan_manual = sm)
#' @references
#' \insertRef{Banerjee:2015}{ssnet}
#'
#' \insertRef{Friedman:2007}{ssnet}
#'
#' \insertRef{Friedman:2010}{ssnet}
#'
#' \insertRef{Morris:2017}{ssnet}
#'
#' \insertRef{Morris:2019}{ssnet}
#'
#' \insertRef{Rockova+George:2018}{ssnet}
#'
#' \insertRef{Tang:2017}{ssnet}
#'
#' @export
#' @inheritParams ssnet_fit
ssnet <- function (x, y, family = c("gaussian", "binomial",
"poisson", "cox"), offset = NULL, epsilon = 1e-04, alpha = 0.5,
maxit = 50, init = NULL, group = NULL, ss = c(0.04, 0.5),
Warning = FALSE, verbose = FALSE, iar.prior = FALSE,
opt.algorithm = "LBFGS", iar.data = NULL, p.bound = c(0.01, 0.99),
tau.prior = "none", stan_manual = NULL, stan_local = FALSE,
plot.pj = FALSE, im.res = NULL)
{
  start.time <- Sys.time()
  call <- match.call()
  x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnames(x) <- paste("x", 1:ncol(x), sep = "")
  nobs <- nrow(x)
  if (NROW(y) != nobs)
    stop("nobs of 'x' and 'y' are different")
  inc <- apply(cbind(y, x), 1, function(z) !any(is.na(z)))
  if (!is.null(offset)) {
    if (length(offset) != nobs)
      stop("nobs of 'x' and 'offset' are different")
    inc <- apply(cbind(y, x, offset), 1, function(z) !any(is.na(z)))
  }
  y <- y[inc]
  x <- x[inc, ]
  offset <- offset[inc]
  family <- family[1]
  if (family == "cox")
    if (!survival::is.Surv(y))
      stop("'y' should be a 'Surv' object")
  if (family == "gaussian")
    y <- (y - mean(y))/sd(y)
  if (!is.null(init) & length(init) != ncol(x))
    stop("give an initial value to each coefficient (not intercept)")
  # create adjacency matrix
  if (iar.prior == TRUE & is.null(iar.data) == TRUE) {
    if (is.null(im.res) == TRUE) {
      stop("Require image dimensions, im.res, to create adjacency matrix. \n")
    } else {
      adjmat <- sim2Dpredictr::binary_adjacency(im.res = im.res, type = "sparse")
      iar.data = mungeCARdata4stan(adjmat$nb.index,
                                   table(adjmat$location.index))
    }
  }
  f <- ssnet_fit(x = x, y = y, family = family, offset = offset, alpha = alpha,
                 epsilon = epsilon, maxit = maxit, init = init, group = group,
                 ss = ss, Warning = Warning, iar.data = iar.data,
                 opt.algorithm = opt.algorithm, p.bound = p.bound,
                 iar.prior = iar.prior, tau.prior = tau.prior,
                 stan_manual = stan_manual, stan_local = stan_local,
                 plot.pj = plot.pj, im.res = im.res)
  f$call <- call
  if (family == "cox")
    class(f) <- c(class(f), "bmlasso", "COXPH")
  else class(f) <- c(class(f), "bmlasso", "GLM")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"),
                   3)
  if (verbose) {
    cat("EM Coordinate Decent Iterations:", f$iter,
        "\n")
    cat("Computational time:", minutes, "minutes \n")
  }
  return(f)
}
