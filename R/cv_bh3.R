#' Cross Validation for Bayesian Models or Elastic Net
#'
#' Several alterations to \code{cv.bh()} were necessary to ensure that \code{update()} works in
#' the functions \code{compare_ssnet()}. Many arguments and functionality are the same as \code{cv.bh()}.
#' See \code{\link[BhGLM]{cv.bh}} for details. An addition in this version is also that for binary
#' outcomes classification and observed accuracy, sensitivity, specificity, and positive and negative
#' predictive values can be output as well as the orginally included measures.
#'
#' @inheritParams BhGLM::cv.bh
#' @inheritParams measure_glm_raw
#' @note The package \code{pROC} will not calculate the AUC when a fold does does not have
#' at least one observation of each level. This can largely be avoided by selecting the number
#' of folds so that such circumstances are rare. When such does occur, the current result is
#' to assign AUC <- NA.
#' @export
cv.bh3 <- function (object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE,
                    classify = FALSE, classify.rule = 0.5)
{
  start.time <- Sys.time()
  out <- cv.bh.lasso2(object = object, nfolds = nfolds,
                      foldid = foldid, ncv = ncv, verbose = verbose,
                      classify = classify, classify.rule = classify.rule)
  stop.time <- Sys.time()
  Time <- round(difftime(stop.time, start.time, units = "min"),
                3)
  if (verbose)
    cat("\n Cross-validation time:", Time, "minutes \n")
  out
}
cv.bh.lasso2 <- function (object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE,
                          classify = FALSE, classify.rule = 0.5)
{
  family <- object$family
  x.obj <- object$x
  y.obj <- object$y
  ss <- object$ss
  group <- object$group
  # n <- NROW(y.obj)
  n <- length(y.obj)
  offset <- object$offset
  init <- object$coefficients
  init <- init[!names(init) %in% "(Intercept)"]
  epsilon <- object$epsilon
  # add for elastic net
  alpha <- object$alpha

  # additional arguments for IAR
  iar.data <- object$iar.data
  p.bound <- object$p.bound
  stan_manual <- object$stan_manual
  stan_local <- object$stan_local
  opt.algorithm <- object$opt.algorithm

  fol <- generate.foldid(nobs = n, nfolds = nfolds, foldid = foldid,
                         ncv = ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- lp0 <- y.fitted0 <- NULL
  j <- 0
  if (verbose)
    cat("Fitting", "ncv*nfolds =", ncv * nfolds,
        "models: \n")
  for (k in 1:ncv) {
    y.fitted <- lp <- rep(NA, n)
    deviance <- NULL
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      if (any(class(object) %in% "glmNet"))
        fit <- stats::update(object, x = x.obj[-omit, ], y = y.obj[-omit],
                      offset = offset[-omit], lambda = object$lambda,
                      alpha = object$alpha,
                      verbose = FALSE)
      if (any(class(object) %in% "bmlasso"))
        fit <- stats::update(object, x = x.obj[-omit, ], y = y.obj[-omit],
                      alpha = object$alpha,
                      offset = offset[-omit], init = init, verbose = FALSE)
      if (is.null(fit$offset))
        fit$offset <- FALSE
      else fit$offset <- TRUE
      xx <- x.obj[omit, , drop = FALSE]
      off <- offset[omit]
      lp[omit] <- as.vector(predict(fit, newx = xx, newoffset = off))
      if (any(class(object) %in% "GLM")) {
        y.fitted[omit] <- as.vector(predict(fit, newx = xx,
                                            type = "response", newoffset = off))
        dd <- suppressWarnings(measure_glm_raw(
                                           y.obj[omit],
                                           y.fitted[omit], family = family, dispersion = fit$dispersion,
                                           classify = classify, classify.rule = classify.rule))
        deviance <- c(deviance, dd["deviance"])
      }
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    if (any(class(object) %in% "GLM")) {
      measures <- measure_glm_raw(y.obj, y.fitted, family = family,
                                  classify = classify, classify.rule = classify.rule)
      measures["deviance"] <- sum(deviance)
      y.fitted0 <- cbind(y.fitted0, y.fitted)
    }
    if (any(class(object) %in% "COXPH"))
      measures <- BhGLM::measure.cox(y.obj, lp)
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
  }
  out <- list()
  if (nrow(measures0) == 1)
    out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE),
                          apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  # out$measures <- round(out$measures, digits = 3)
  out$y.obs <- y.obj
  out$lp <- lp0
  if (any(class(object) %in% "GLM"))
    out$y.fitted <- y.fitted0
  out$foldid <- foldid
  out
}
