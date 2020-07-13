#' Fit the Spike-and-Slab Elastic Net GLM with Intrinsic Autoregressive Prior
#'
#' @importFrom rstan optimizing
#' @importFrom stats coef deviance predict sd
#' @importFrom survival Surv is.Surv
#' @importFrom glmnet glmnet
#' @importFrom BhGLM bglm
#' @param x Design, or input, matrix, of dimension nobs x nvars; each row is an observation vector. It is recommended that
#' \code{x} have user-defined column names for ease of identifying variables. If missing \code{colnames} are internally
#' assigned \code{x1}, \code{x2}, ... and so forth.
#' @param y Scalar response variable. Quantitative for \code{family = "gaussian"}, or \code{family = "poisson"}
#' (non-negative counts). For \code{family = "gaussian"}, \code{y} is always standardized. For \code{family = "binomial"},
#' y should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated
#' as the target class; for a factor, the last level in alphabetical order is the target class). For \code{family="cox"}, \code{y}
#' should be a two-column matrix with columns named \code{'time'} and \code{'status'}. The latter is a binary variable, with \code{'1'}
#' indicating death, and \code{'0'} indicating right censored. The function \code{Surv()} in package survival produces such a matrix.
#' @param family Response type (see above).
#' @param alpha A scalar value between 0 and 1 determining the compromise between the Ridge and Lasso models. When
#' \code{alpha = 1} reduces to the Lasso, and when \code{alpha = 0} reduces to Ridge.
#' @param iar.prior Logical. When \code{TRUE}, imposes intrinsic autoregressive prior on logit of the probabilities of
#' inclusion. When \code{FALSE}, treats probabilties of inclusion as unstructured.
#' @param offset A vector of length \code{nobs} that is included in the linear predictor.
#' @param epsilon A positive convergence tolerance; the iterations converge when \eqn{|dev - dev_old|/(|dev| + 0.1) < e}.
#' @param maxit An integer giving the maximal number of EM iterations.
#' @param init A vector of initial values for all coefficients (not for intercept). If not given, it will be internally produced.
#' @param group A numeric vector, or an integer, or a list indicating the groups of predictors.
#' If \code{group = NULL}, all the predictors form a single group. If \code{group = K}, the
#' predictors are evenly divided into groups each with K predictors. If group is a numberic vector,
#' it defines groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2: \code{(group[2]+1):group[3]},
#' Group 3: \code{(group[3]+1):group[4]}, ... If group is a list of variable names, \code{group[[k]]}
#' includes variables in the k-th group. The mixture double-exponential prior is only used for grouped
#' predictors. For ungrouped predictors, the prior is double-exponential with scale \code{ss[2]} and mean 0.
#' @param ss A vector of two positive scale values for the spike-and-slab mixture double-exponential prior, allowing for different
#' scales for different predictors, leading to different amount of shrinkage. Smaller scale values give stronger shrinkage. While
#' the smaller of the two input values will be treated as the spike scale, it is recommended to specify the spike scale as the
#' first element of the vector.
#' @param Warning Logical. If \code{TRUE}, shows the error messages of not convergence and identifiability.
#' @param adjmat A data.frame or matrix containing a "sparse" representation of the neighbor relationships. The first
#' column should contain a numerical index for a given location. Each index will be repeated in this column for
#' every neighbor it has. The indices for the location's neighbors are then specified in the second column. An optional
#' third column specifies weights. If no third column is specified, then equal weights are assumed.
#' @param iar.data A list of output from \code{\link{mungeCARdata4stan}} that contains the necessary
#' inputs for the IAR prior. When unspecified, this is built internally assuming that neighbors are those
#' variables directly above, below, left, and  right of a given variable location. \code{im.res} must be specified
#' when allowing this argument to be built internally. It is not recommended to use this argument directly, even
#' when specifying a more complicated neighborhood stucture; this can be specified with the \code{adjmat} argument,
#' and then internally converted to the correct format.
#' @param opt.algorithm One of \code{c("LBFGS", "BFGS", "Newton")}. This argument determines which argument
#' is used to optimize the term in the EM algorithm that estimates the probabilities of inclusion for
#' each parameter. Optimization is performed by \code{optimizing}.
#' @param tau.prior One of \code{c("none", "manual", "cauchy")}. This argument determines the precision
#' parameter in the Conditional Autoregressive model for the (logit of) prior inclusion probabilities.
#' When \code{"none"}, the precision is set to 1; when "manual", the precision is manually entered by the
#' use; when \code{"cauchy"}, the inverse precision is assumed to follow a Cauchy distribution with mean 0 and
#' scale 2.5.
#' @param plot.pj When \code{TRUE}, prints a series of 2D graphs of the prior probabilities of inclusion
#' at each step of the algorithm. This should NOT be used for 3D data.
#' @param im.res A 2-element vector where the first argument is the number of "rows" and the second argument
#' is the number of "columns" in each subject's "image". Default is \code{NULL}.
#' @param p.bound A vector defining the lower and upper boundaries for the probabilities of inclusion
#' in the model, respectively. Defaults to \code{c(0.01, 0.99)}.
#' @param stan_manual A \code{stan_model} that is manually specified. Especially when fitting multiple models in
#' succession, specifying the the \code{stan} model outside this "loop" may avoid errors.
#' @param stan_local Logical. Defaults to \code{FALSE}, but when \code{TRUE}, uses locally stored \code{stan} files.
#' This option will eventually be removed once the package is more stable.
#' @return The fitted model for the spike-and-slab elastic net. An object of class \code{c("elnet", "glmnet"}.
#' @note Currently, the \code{ssnet()} \code{im.res} can only handle 2D data. Future versions may allow
#' images to be 3D. However, the function will work given any appropriately specified neighborhood matrix,
#' whatever the original dimension.
#' @examples
#' cn <- c()
#' for (i in 1:100) cn[i] <- paste0("x", i)
#' bml_model <- ssnet_fit(x = matrix(rnorm(10000), nrow = 100, ncol = 100,
#'                              dimnames = list(1:100, cn)), family = "gaussian",
#'                              y = rnorm(100))
ssnet_fit <- function (x, y, family = c("gaussian", "binomial", "poisson", "cox"),
                       offset = NULL, epsilon = 1e-04, alpha = 0.5,
                       maxit = 50, init = rep(0, ncol(x)), group = NULL,
                       ss = c(0.04,0.5), Warning = FALSE,
                       iar.prior = FALSE, adjmat = NULL,  iar.data = NULL,
                       tau.prior = "none", stan_manual = NULL, stan_local = FALSE,
                       opt.algorithm = "LBFGS", p.bound = c(0.01, 0.99),
                       plot.pj = FALSE, im.res = NULL)
{
  if (plot.pj == TRUE & !requireNamespace("sim2Dpredictr") == FALSE) {
    stop("Cannot plot p_j without package sim2Dpredictr. \n")
  }

  ss <- sort(ss)
  prior.scale <- ss[length(ss)]
  if (family == "cox")
    intercept <- FALSE
  else intercept <- TRUE
  x0 <- x
  if (intercept)
    x0 <- cbind(1, x)
  d <- prepare(x = x0, intercept = intercept, prior.mean = 0,
               prior.sd = 1, prior.scale = prior.scale, prior.df = 1,
               group = group)
  x <- d$x
  prior.scale <- d$prior.scale
  group <- d$group
  group.vars <- d$group.vars
  ungroup.vars <- d$ungroup.vars
  prior.scale <- prior.scale/autoscale(x, min.x.sd = 1e-04)

  if (intercept) {
    x <- x[, -1]
    prior.scale <- prior.scale[-1]
  }

  gvars <- unlist(group.vars)
  theta <- p <- rep(0.5, length(gvars))
  names(theta) <- names(p) <- gvars
  if (is.null(init)) {
    for (k in 1:5) {
      ps <- ss[1] + (k - 1) * 0.01
      if (family == "cox")
        ps <- min(ss[1] + (k - 1) * 0.01, 0.08)
      f <- glmnet::glmnet(x = x, y = y, family = family, offset = offset,
                          alpha = 0.95, lambda = 1/(nrow(x) * ps), standardize = TRUE)
      b <- as.numeric(f$beta)
      if (any(b != 0))
        break
    }
  }
  else b <- as.numeric(init)
  names(b) <- colnames(x)
  b <- ifelse(b == 0, 0.001, b)
  init <- b
  devold <- 0
  conv <- FALSE
  for (iter in 1:maxit) {
    out <- update_scale_p(b0 = b[gvars], ss = ss, theta = theta)
    prior.scale[gvars] <- out[[1]]
    p <- out[[2]]
    if (plot.pj == TRUE) {
      if (is.null(im.res) == TRUE) {
        stop("Require image dimensions, im.res, to produce image plot. \n")
      }
      cat("The number of non-zero parameters is ", length(b[b != 0]), "\n")
      sim2Dpredictr::inf_2D_image(B = p, im.res = im.res,
                                  binarize.B = FALSE,
                                  B.incl.B0 = FALSE)
    }
    if (iar.prior == FALSE) {
      theta <- update_ptheta_group2(group.vars = group.vars,
                                    p = p, p.bound = p.bound)
    } else {
      theta <- max_q2_iar(iar.data = iar.data, p = p,
                          opt.algorithm = opt.algorithm,
                          tau.prior = tau.prior,
                          stan_manual = stan_manual,
                          stan_local = stan_local,
                          p.bound = p.bound)
    }
    Pf <- 1/(prior.scale + 1e-10)
    f <- glmnet::glmnet(x = x, y = y, family = family, offset = offset,
                        alpha = alpha, penalty.factor = Pf,
                        lambda = sum(Pf)/(nrow(x) * ncol(x)),
                        standardize = FALSE)
    b <- as.numeric(f$beta)
    names(b) <- colnames(x)
    dev <- deviance(f)
    if (abs(dev - devold)/(0.1 + abs(dev)) < epsilon & iter >
        5) {
      conv <- TRUE
      break
    }
    else devold <- dev
  }
  if (Warning & !conv)
    warning("algorithm did not converge", call. = FALSE)
  f$x <- x
  f$y <- y
  f$family <- family
  f$ss <- ss
  f$coefficients <- as.numeric(coef(f))
  names(f$coefficients) <- rownames(coef(f))
  f$linear.predictors <- predict(f, newx = x, type = "link",
                                 offset = offset)
  if (family == "gaussian")
    f$dispersion <- BhGLM::bglm(y ~ f$linear.predictors - 1, start = 1,
                         prior = "de", prior.mean = 1, prior.scale = 0,
                         verbose = FALSE)$dispersion
  f$iter <- iter
  f$prior.scale <- prior.scale
  f$penalty.factor <- Pf
  f$group <- group
  f$group.vars <- group.vars
  f$ungroup.vars <- ungroup.vars
  f$p <- p
  f$ptheta <- theta
  f$init <- init
  f$aic <- deviance(f) + 2 * f$df
  f$offset <- offset
  f$iar.data <- iar.data
  f$stan_manual <- stan_manual
  f$stan_local <- stan_local
  f$opt.algorithm <- opt.algorithm
  f$epsilon <- epsilon
  f$alpha <- alpha
  return(f)
}
