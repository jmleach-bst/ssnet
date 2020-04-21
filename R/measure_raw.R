#' Evaluating Fitted Models
#'
#' @note This function is taken directly from \code{measure.glm} in \code{BhGLM}, with the modification that
#' measures are no longer rounded.
#'
#' @param y This is an outcome/response vector.
#' @param y.fitted This predicted (estimated) response values for GLMs or probabilties of
#' response values for ordinal models from a fitted model or cross-validation.
#' @param family A character stating to which family the model belongs.
#' @param dispersion A scalar defining the dispersion parameter from a GLM, or \eqn{theta}
#' for negative binomial.
#' @return A vector.
measure_glm_raw <- function(y, y.fitted, family, dispersion = 1)
{
  if (NROW(y) != NROW(y.fitted))
    stop("y and y.fitted should be of the same length",
         call. = FALSE)
  if (is.null(dispersion))
    dispersion <- 1
  mu <- y.fitted
  if (substr(family, 1, 6) == "NegBin" | substr(family,
                                                1, 17) == "Negative Binomial" | family == "nb")
    family <- "NegBin"
  fam <- c("gaussian", "binomial", "poisson",
           "quasibinomial", "quasipoisson", "NegBin")
  if (!family %in% fam)
    stop("Measures for this family have not been implemented yet")
  if (family == "gaussian")
    logL <- dnorm(y, mean = mu, sd = sqrt(dispersion), log = TRUE)
  if (family == "binomial" | family == "quasibinomial") {
    if (is.factor(y))
      y <- as.numeric(y) - 1
    L <- dbinom(y, size = 1, prob = mu, log = FALSE)
    L <- ifelse(L == 0, 1e-04, L)
    logL <- log(L)
  }
  if (family == "poisson" | family == "quasipoisson")
    logL <- dpois(y, lambda = mu, log = TRUE)
  if (family == "NegBin")
    logL <- dnbinom(y, size = dispersion, mu = mu, log = TRUE)
  logL <- sum(logL, na.rm = TRUE)
  deviance <- -2 * logL
  mse <- mean((y - mu)^2, na.rm = TRUE)
  mae <- mean(abs(y - mu), na.rm = TRUE)
  measures <- list(deviance = deviance, mse = mse, mae = mae)
  if (family == "gaussian") {
    R2 <- (var(y, na.rm = TRUE) - mse)/var(y, na.rm = TRUE)
    measures <- list(deviance = deviance, R2 = R2, mse = mse,
                     mae = mae)
  }
  if (family == "binomial" | family == "quasibinomial") {
    if (!requireNamespace("pROC"))
      install.packages("pROC")
    require(pROC)
    # print(y)
    if (length(unique(y)) != 2) {
      # at present we cannot obtain AUC when all(y) == 1 or all(y) == 0
      AUC <- NA
      # cat("Outcomes for this fold are either all 0 or all 1. Returning auc = NA. \n")
    } else {
      AUC <- suppressMessages(pROC::auc(y, mu))
      AUC <- as.numeric(AUC)
    }
    misclassification <- mean(abs(y - mu) >= 0.5, na.rm = TRUE)
    measures <- list(deviance = deviance, auc = AUC, mse = mse,
                     mae = mae, misclassification = misclassification)
  }
  # round(unlist(measures), digits = 3)
  unlist(measures)
}

#' Obtain Measures of Model Fitness
#'
#' @note This function is taken directly from \code{measure.bh} in \code{BhGLM}, with the modification that
#' measures are no longer rounded.
#'
#' @param object A fitted object..
#' @param new.x A data frame or matrix of new values for variables used in object.
#' For an object from \code{bglm}, \code{bpolr} or \code{bcoxph}, it is data frame;
#' for object from \code{glmNet} or \code{bmlasso}, it is matrix.
#' @param new.y A vector of new response values corresponding to new.x. If \code{new.x}
#' or \code{new.y} are omitted, the fitted linear predictors are used for prediction.
#' @param new.offset A data frame or vector of offset values for new data points.
#' If \code{new.x} includes offset, do not need to set \code{new.offset}.
#' @return A vector.
measure_bh_raw <- function (object, new.x, new.y, new.offset)
{
  if (missing(new.y)) {
    y <- object$y
    if (any(class(object) %in% "polr"))
      y <- object$model[, 1]
  }
  else y <- new.y
  if (!missing(new.x))
    if (NROW(new.x) != NROW(y))
      stop("'new.x' and 'new.y' should be the same length.",
           call. = FALSE)
  if (missing(new.offset))
    new.offset <- NULL
  if (!is.null(new.offset))
    if (NROW(new.offset) != NROW(y))
      stop("'new.offset' and 'new.y' should be the same length.",
           call. = FALSE)
  if (!is.null(new.offset) & !missing(new.x) & any(class(object) %in%
                                                   c("glm", "coxph", "polr")))
    new.x <- data.frame(new.x, new.offset)
  if (any(class(object) %in% "glm")) {
    mu <- predict(object, newdata = new.x, type = "response")
    if (any(class(object) %in% "negbin"))
      object$dispersion <- object$theta
    measures <- measure.glm(y, mu, family = object$family$family,
                            dispersion = object$dispersion)
  }
  if (any(class(object) %in% "glmNet") | any(class(object) %in%
                                             "bmlasso")) {
    family <- object$family
    if (family == "cox")
      if (!is.Surv(y))
        stop("'new.y' must be a Surv object")
    if (is.null(object$offset))
      object$offset <- FALSE
    else object$offset <- TRUE
    if (missing(new.x))
      new.x <- object$x
    newx <- as.matrix(new.x)
    if (family == "cox") {
      lp <- predict(object, newx = newx, newoffset = new.offset)
      lp <- as.vector(lp)
      measures <- measure.cox(y, lp)
    }
    else {
      mu <- predict(object, newx = newx, newoffset = new.offset,
                    type = "response")
      mu <- as.vector(mu)
      measures <- measure_glm_raw(y, mu, family = family, dispersion = object$dispersion)
    }
  }
  if (any(class(object) %in% "coxph")) {
    if (!is.Surv(y))
      stop("'new.y' must be a Surv object")
    lp <- predict(object, newdata = new.x)
    measures <- measure.cox(y, lp)
  }
  if (any(class(object) %in% "polr")) {
    if (!is.factor(y))
      stop("'new.y' must be a factor")
    probs <- predict(object, newdata = new.x, type = "probs")
    measures <- measure.polr(y, probs)
  }
  measures
}
