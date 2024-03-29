#' Extract True/False Positive/Negative Measures After Classification
#'
#' After fitting GLM's with, set a classification  rule, and obtain measures
#' for classification effectiveness.
#'
#' @param num.class A scalar value defining the number of classes. Default
#' value is two.
#' @param classify.rule A vector of thresholds that determine assigned classes.
#' @param beta.hat A vector of parameter estimates from some model fit.
#' @param x A (numeric) design matrix whose number of columns is equal to
#' \code{length(beta.hat) - 1} and whose number of rows is equal to
#' \code{length(y)}. Strictly speaking, \code{x} may be a data frame, but
#' categorical variables must have already been transformed into alternative
#' coding schemes; e.g., dummy coding.
#' @param make.x0 Logical. When \code{TRUE}, creates a column of ones for the
#' intercept. Default is \code{TRUE}.
#' @param y A numeric vector of outcome measurements. The length of this
#' vector must equal \code{nrow(x)}.
#' @param family One of \code{"gaussian", "gamma", "poisson", "binomial"},
#' which determines the distribution from which \code{y} is assumed to be
#' drawn.
#' @param gamma.link Defines the link function in gamma GLM's. One of
#' \code{c("inverse", "log", "identity")}. The default is \code{"log"}.
#' @param extra.measures Logical. When \code{TRUE}, also returns mean squared
#' error (mse) and mean absolute error (mae), and when
#' \code{family = "binomial"}, returns AUC and misclassification. Default is
#' \code{FALSE}.
#' @return A data frame with a single row containing columns for sensitivity,
#' specificity, ppv, npv, and accuracy, where ppv is positive predictive value
#' and npv is negative predictive value. When there are 3 or more classes,
#' only accuracy, and possibly mse and mae, are returned.
#' @details Note that \code{y} and \code{x} may be the same outcomes/design
#' matrix used to generate \code{beta.hat}, but they may also, and perhaps
#' more appropriately, be a held-out/independent data set on which to test the
#' classification performance of whatever model was used to obtain
#' \code{beta.hat}.
#' @examples
#' set.seed(72874)
#'
#' # sample size
#' n <- 10
#'
#' # true betas
#' b <- c(1, 0.1, 0.5, 0, -0.1, -0.5)
#'
#' # observed predictors
#' x <- matrix(rnorm(n * (length(b) - 1)), nrow = n, ncol = length(b) - 1)
#' x0 <- rep(1, n)
#'
#' # observed outcomes
#' y <- cbind(x0, x) %*% b + rnorm(n, 0, 2.5)
#' y <- ifelse(y > 0, 1, 0)
#'
#' # fitted model
#' yx.df <- data.frame(y, x)
#' beta.hat <- glm(y ~ ., data = yx.df , family = "binomial")$coefficients
#'
#' # evaluate classification
#' eval_classify(classify.rule = 0.5, beta.hat = beta.hat,
#'               x = x, y = y, family = "binomial",
#'               extra.measures = TRUE)
#' @export
eval_classify <- function(
    num.class = 2,
    classify.rule,
    beta.hat,
    x,
    y,
    family,
    make.x0 = TRUE,
    gamma.link = "log",
    extra.measures = FALSE
){
  # checks
  if (num.class < 2) {
    stop("Require at least 2 classes.")
  }

  if (is.numeric(beta.hat) == FALSE) {
    if (is.numeric(t(beta.hat)) == TRUE) {
      beta.hat <- t(beta.hat)
    } else {
      stop("Parameter vector beta.hat must contain only numeric values.")
    }
  }

  if (is.numeric(y) == FALSE) {
    if (is.numeric(t(y)) == TRUE) {
      y <- t(y)
    } else {
      stop("y must contain only numeric values.")
    }
  }

  if (is.numeric(x) == FALSE) {
    if (is.numeric(as.matrix(x)) == TRUE) {
      x <- as.matrix(x)
    } else {
      stop("x must contain only numeric values.")
    }
  }

  # create vector of ones for intercept when necessary
  if (make.x0 == TRUE) {
    x0 <- rep(1, nrow(x))
    x <- cbind(x0, x)
  }

  # calculate predicted values/expectations
  if (family == "gaussian" | (family == "gamma" & gamma.link == "identity")) {
    y.hat <- x %*% beta.hat
  }

  if (family == "binomial") {
    y.hat <- 1 / (1 + 1/(exp(x %*% beta.hat)))
  }

  if (family == "poisson" | (family == "gamma" & gamma.link == "log")) {
    y.hat <- exp(x %*% beta.hat)
  }

  if (family == "gamma" & gamma.link == "inverse") {
    y.hat <- 1 / (x %*% beta.hat)
  }

  # obtain classifications
  classy <- function(outcome, classify.rule) {
    classify <- rep(0, length(outcome))
    n.cat <- length(classify.rule) + 1
    for (i in 1:n.cat) {
      if (i == 1) {
        classify[outcome > classify.rule[1]] <- 1
      }
      if (i > 1 & i < n.cat) {
        classify[(classify.rule[i] < outcome) &
                   (outcome < classify.rule[i+1])] <- i
      }
      if (i == n.cat) {
        classify[outcome > classify.rule[i - 1]] <- n.cat - 1
      }
    }
    return(classify)
  }

  predicted.class <- classy(
    outcome = y.hat,
    classify.rule = classify.rule
  )
  observed.class <- classy(
    outcome = y,
    classify.rule = classify.rule
  )

  # obtain classification performance measures
  accuracy <- length(which(predicted.class == observed.class)) / length(y)

  if (num.class == 2) {

    # calculate true/false positive/negative
    fp <- rep(0, length(y))
    tp <- rep(0, length(y))
    fn <- rep(0, length(y))
    tn <- rep(0, length(y))

    fp[predicted.class == 1 & observed.class == 0] <- 1
    tp[predicted.class == 1 & observed.class == 1] <- 1
    fn[predicted.class == 0 & observed.class == 1] <- 1
    tn[predicted.class == 0 & observed.class == 0] <- 1

    # number of observed positives/negatives
    op <- sum(observed.class)
    on <- length(observed.class) - op

    # number of predicted positives/negatives
    pp <- sum(predicted.class)
    pn <- length(predicted.class) - pp

    sensitivity <- sum(tp) / op
    specificity <- sum(tn) / on
    ppv <- sum(tp) / pp
    npv <- sum(tn) / pn

    predm <- data.frame(
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      ppv = ppv,
      npv = npv
    )
  } else {
    predm <- data.frame(
      accuracy == accuracy
    )
  }

  if (extra.measures == TRUE) {
    # mean square error
    mse <- mean((y - y.hat) ^ 2)

    # mean absolute error
    mae <- mean(abs(y - y.hat))

    m.df <- data.frame(
      mse = mse,
      mae = mae
    )

    if (family == "binomial") {
      # auc
      yhatn <- as.numeric(y.hat)
      yn <- as.numeric(y)
      AUC <- suppressMessages(pROC::auc(yn, yhatn))
      auc <- as.numeric(AUC)

      # misclassification
      ae <- abs(y - y.hat)
      aei <- rep(0, length(y))
      aei[ae > 0.5] <- 1
      misclassification <- mean(aei)
    }

    predm <- cbind(
      auc,
      m.df,
      misclassification,
      predm
    )
  }

  return(predm)
}















