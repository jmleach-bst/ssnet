#' Bayesian Spike-and-Slab Elastic Net with Spatial Structure
#'
#' Fits generalized linear models with spike-and-slab priors whose (logit of) probability of inclusion in the model
#' is either assigned an intrinsic autoregression as a prior to incorporate spatial information or is unstructured.
#' The model is fit using an EM algorithm where the E-step is fit by \code{glmnet()} and the M-step is fit using the
#' \code{stan} function \code{optimizing} (when IAR prior is employed).
#'
#' @inheritParams ssnet_fit
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
#' @note While the type.multinomial option is included, it is only valid for traditional elastic net models.
#' Thus far we have only extended the spike-and-slab models for grouped selection.
#' @examples
#' library(sim2Dpredictr)
#' set.seed(223)
#'
#' ## sample size
#' n <- 30
#' ## image dims
#' nr <- 4
#' nc <- 4
#'
#' ## generate data
#' cn <- paste0("x", 1:(nr * nc))
#' tb <- rbinom(nr * nc, 1, 0.2)
#' tx <- matrix(rnorm(n * nr * nc), nrow = n, ncol = nr * nc,
#'              dimnames = list(1:n, cn))
#' ty <- tx %*% tb + rnorm(n)
#'
#' ## build adjacency matrix
#' adjmat <- proximity_builder(im.res = c(nr, nc), type = "sparse")
#' ## stan model information
#' model_info <- mungeCARdata4stan(adjmat$nb.index,
#'                                 table(adjmat$location.index))
#'
#' ## fit model
#' ex_model <- ssnet(x = tx, y = ty, alpha = 0.5,
#'                   iar.prior = TRUE, iar.data = model_info,
#'                   family = "gaussian")
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
#' \insertRef{Tang:2018}{ssnet}
#'
#' @export
ssnet <- function (x, y, family = c("gaussian", "binomial", "multinomial", "poisson", "cox"),
                   offset = NULL, epsilon = 1e-04, alpha = 0.95, type.multinomial = "grouped",
                   maxit = 50, init = rep(0, ncol(x)), init.theta = 0.5,
                   ss = c(0.04,0.5), Warning = FALSE, group = NULL,
                   iar.prior = FALSE, adjmat = NULL,  iar.data = NULL,
                   tau.prior = "none", tau.manual = NULL, stan_manual = NULL,
                   opt.algorithm = "LBFGS", p.bound = c(0.01, 0.99),
                   plot.pj = FALSE, im.res = NULL, verbose = FALSE, print.iter = FALSE)
{
  start.time <- Sys.time()
  call <- match.call()

  ##############
  # check data #
  ##############

  # picked a family?
  if (length(family) != 1){
    stop("user must select a family: gaussian, binomial, multinomial, poisson, or cox")
  }

  # check matrix
  if (is.matrix(x) == FALSE) {
    stop("x should be a matrix")
  }

  # assign variable names if necessary
  if (is.null(colnames(x)) == TRUE) {
    colnames(x) <- paste0("x", 1:ncol(x))
  }

  # ensure number obs in x and y are the same
  if (is.null(length(y)) == FALSE) {
    if (length(y) != nrow(x)) {
      stop("length of y should equal number of rows in x")
    }
  } else {
    if (nrow(y) != nrow(x)) {
      stop("length of y should equal number of rows in x")
    }
  }

  # check that y & family makes sense
  if (is.numeric(y) == FALSE & (family == "gaussian") | (family == "poisson")) {
    stop("gaussian family requires numeric y")
  }

  # check EN parameter
  if (alpha > 1 | alpha < 0) {
    stop("alpha must be 0 or greater and cannot exceed 1.")
  }

  # check scale parameters
  if (any(ss <= 0)) {
    stop("scale values must exceed 0.")
  }

  if (length(ss) != 2) {
    stop("ss should contain only 2 elements, 1 spike and 1 slab scale.")
  }

  nobs <- nrow(x)

  # remove observations with missingness
  inc <- apply(cbind(y, x), 1, function(z) !any(is.na(z)))
  if (!is.null(offset)) {
    if (length(offset) != nobs) {
      stop("nobs of 'x' and 'offset' are different")
    }
        inc <- apply(cbind(y, x, offset), 1, function(z) !any(is.na(z)))
  }
  y <- y[inc]
  x <- x[inc, ]
  offset <- offset[inc]

  # cox requires y to be a survival object
  if (family == "cox")
    if (!survival::is.Surv(y))
      stop("'y' should be a 'Surv' object")

  # standardize y for continuous data
  if (family == "gaussian") {
    y <- (y - mean(y)) / sd(y)
  }

  # when specifying initial values, ensure there are the right number
  if (!is.null(init) & length(init) != ncol(x))
    stop("must specify initial value to each coefficient except intercept")

  # format neighborhood information
  if (is.null(adjmat) == FALSE & is.null(iar.data) == FALSE) {
    warning("argument adjmat is not used when argument iar.data is specified.")
  }

  # format IAR data
  if (iar.prior == TRUE & is.null(iar.data) == TRUE) {
    iar.data <- format_iar(adjmat = adjmat, im.res = im.res, x = x)
  }

  f <- ssnet_fit(x = x, y = y, family = family, offset = offset, alpha = alpha,
                 epsilon = epsilon, maxit = maxit, init = init, init.theta = init.theta,
                 group = group, ss = ss, Warning = Warning, iar.data = iar.data,
                 opt.algorithm = opt.algorithm, p.bound = p.bound,
                 iar.prior = iar.prior, tau.prior = tau.prior,
                 stan_manual = stan_manual,
                 plot.pj = plot.pj, im.res = im.res,
                 print.iter = print.iter)
  f$call <- call
  if (family == "cox") {
    class(f) <- c(class(f), "bmlasso", "COXPH")
  } else {
    class(f) <- c(class(f), "bmlasso", "GLM")
  }
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
