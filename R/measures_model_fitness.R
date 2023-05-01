#' Obtain Measures of Model Fitness
#'
#' TBD
#'
#' @param y A vector containing responses/outcome values..
#' @param y.fitted A vector containing predicted outcome/response values
#' obtain from some model.
#' @param inverse.link.y Logical. When \code{TRUE}, \code{y.fitted} is assumed
#' to have had the inverse link function applied to the values. When
#' \code{FALSE}, \code{y.fitted} is assumed to be the linear predictor
#' \eqn{XB}, and the inverse link function is internally applied to
#' \code{y.fitted}.
#' @param classify Logical. When \code{TRUE} and \code{family = "binomial"}
#' applies a classification rule given by the argument \code{classify.rule},
#' and outputs accuracy, sensitivity, specificity, positive predictive value
#' (ppv), and negative predictive value (npv).
#' @param classify.rule A value between 0 and 1. For a given predicted value
#' from a logistic regression, if the value is above \code{classify.rule},
#' then the predicted class is 1; otherwise the predicted class is 0. The
#' default is 0.5.
#' @param dispersion A scalar defining the dispersion parameter from a GLM,
#' or \eqn{theta} for negative binomial.
#' @inheritParams ssnet
#' @return A data frame consisting of a single row and a column for each model
#' fitness measure.
#' @note This function is set to replace \code{measure_glm_raw} and
#' \code{measure_bh_raw}.
#' @details When the family is treated as Gaussian, returns deviance, R2, mean
#' squared error (MSE), and mean absolute error (MAE). Additionally, when the
#' outcome is binary, returns misclassification, and if \code{classify = TRUE},
#' then returns accuracy, sensitivity, specificity, positive predictive value
#' (PPV), negative predictive value (NPV), Matthews correlation coefficient
#' (MCC), and F1 score.
#' @examples
#'
#' ## binary data
#' yb <- rbinom(20, size = 1, prob = 0.5)
#' yb.f <- yb
#' yb.f[c(2, 5, 6, 7, 13, 17)] <- abs(1 - yb[c(2, 5, 6, 7, 13, 17)])
#'
#' measures_model_fitness(
#'   y = yb,
#'   y.fitted = yb.f,
#'   family = "binomial",
#'   classify = TRUE
#'   )
#'
#' ## gaussian data
#' yg <- rnorm(20)
#' yg.f <- yg + rnorm(20, 0, 2/3)
#'
#' measures_model_fitness(
#'   y = yg,
#'   y.fitted = yg.f,
#'   family = "gaussian",
#'   dispersion = 1
#'  )
#'
#' @export
measures_model_fitness <- function(
  y,
  y.fitted,
  family,
  dispersion = NULL,
  inverse.link.y = TRUE,
  classify = FALSE,
  classify.rule = 0.5
){
  if (length(y) != length(y.fitted)) {
    stop("y and y.fitted should be of equal length")
  }

  if (!(family %in% c("gaussian", "binomial", "poisson", "NegBin"))) {
    stop("family must be one of gaussian, binomial, NegBin, or poisson.")
  }

  if (is.null(dispersion)) {
    if (family %in% c("binomial", "poisson")) {
      dispersion <- 1
    } else {
      warning(
        "No dispersion parameter provided. \n
         Cannot estimate loglikelihood or deviance."
        )
    }
  }

  # calculate log likelihood and deviance
  if (family == "gaussian") {
    if (is.null(dispersion)) {
      logL <- NULL
    } else {
      logL <- stats::dnorm(
        y,
        mean = y.fitted,
        sd = sqrt(dispersion),
        log = TRUE
      )
    }
  }
  if (family == "binomial") {
    if (is.factor(y)) {
      y <- as.numeric(y) - 1
    }
    if (inverse.link.y == FALSE) {
      y.fitted <- exp(y.fitted) / (1 + exp(y.fitted))
    }
    L <- stats::dbinom(
      y,
      size = 1,
      prob = y.fitted,
      log = FALSE
    )
    L <- ifelse(L == 0, 1e-04, L)
    logL <- log(L)
  }
  if (family == "poisson") {
    if (inverse.link.y == FALSE) {
      y.fitted <- exp(y.fitted)
    }
    logL <- stats::dpois(
      y,
      lambda = y.fitted,
      log = TRUE
    )
  }

  if (family == "NegBin") {
    if (inverse.link.y == FALSE) {
      y.fitted <- exp(y.fitted)
    }
    logL <- stats::dbinom(
      y,
      size = dispersion,
      mu = y.fitted,
      log = TRUE
    )
  }

  logL <- sum(logL, na.rm = TRUE)
  deviance <- -2 * logL

  # Mean Squared & Absolute Errors
  mse <- mean((y - y.fitted) ^ 2, na.rm = TRUE)
  mae <- mean(abs(y - y.fitted), na.rm = TRUE)

  # if (!(family %in% c("gaussian", "binomial"))) {
  #   measures <- data.frame(
  #     deviance = deviance,
  #     mse = mse,
  #     mae = mae
  #   )
  # }

  # Additional measures for Gaussian data
  if (family == "gaussian") {
    R2 <- (stats::var(y, na.rm = TRUE) - mse) / stats::var(y, na.rm = TRUE)
    measures <- data.frame(
      deviance = deviance,
      R2 = R2,
      mse = mse,
      mae = mae
    )
  }

  # Additional measures of Binomial data
  if (family == "binomial") {
    # if (!requireNamespace("pROC"))
    #   install.packages("pROC")
    # require(pROC)
    # print(y)
    if (length(unique(y)) != 2) {
      # at present we cannot obtain AUC when all(y) == 1 or all(y) == 0
      AUC <- NA
      # cat("Outcomes for this fold are either all 0 or all 1. Returning auc = NA. \n")
    } else {
      y.roc <- as.numeric(y)
      y.fitted.roc <- as.numeric(y.fitted)
      AUC <- suppressMessages(pROC::auc(y.roc, y.fitted.roc))
      AUC <- as.numeric(AUC)
    }
    misclassification <- mean(abs(y - y.fitted) >= 0.5, na.rm = TRUE)
    measures <- data.frame(
      deviance = deviance,
      auc = AUC,
      mse = mse,
      mae = mae,
      misclassification = misclassification
    )

    if (classify == TRUE) {
      predicted.class <- rep(0, length(y))
      predicted.class[y.fitted > classify.rule] <- 1
      observed.class <- y

      # Matthew's Correlation Coefficient
      mcc <- function(tp, tn, fp, fn) {
        mcc.num <- (tp * tn) - (fp * fn)
        mcc.den <- (tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)
        return(mcc.num / sqrt(mcc.den))
      }

      # f1 score
      f1.score <- function(tp, tn, fp, fn) {
        f1.num <- 2 * tp
        f1.den <- 2 * tp + fn + fp
        return(f1.num / f1.den)
      }

      # obtain classification performance measures
      accuracy <- length(which(predicted.class == observed.class)) / length(y)

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
      mcc.obs <- mcc(
        tp = sum(tp),
        tn = sum(tn),
        fp = sum(fp),
        fn = sum(fn))
      f1.obs <- f1.score(
        tp = sum(tp),
        tn = sum(tn),
        fp = sum(fp),
        fn = sum(fn)
      )

      predm <- data.frame(
        accuracy = accuracy,
        sensitivity = sensitivity,
        specificity = specificity,
        ppv = ppv,
        npv = npv,
        mcc = mcc.obs,
        f1 = f1.obs
      )

      measures <- cbind(measures, predm)
    }
  }
  return(measures)
}
