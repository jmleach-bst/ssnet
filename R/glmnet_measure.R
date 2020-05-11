#' Extract Measures of Model Fit from GLMNET Fit
#'
#' Extracts measures of model fitness from a \code{glmnet} object.
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @param cv.model A model fit by \code{glmnet} function \code{cv.glmnet()}.
#' @param lambda.criteria Determines the model selection criteria. When \code{"lambda.min"} the
#' final model is selected based on the penalty that minimizes the measure given in \code{type.measure}.
#' When \code{"lambda.1se"} the final model is selected based on the smallest value of lambda that
#' is within one standard error of the minimal measure given in \code{type.measure}.
#' @inheritParams glmnet::glmnet
#' @inheritParams glmnet::cv.glmnet
#' @return A data frame containing measures of model fitness.
#' @note This function is primarily used within \code{\link[ssnet]{compare_ssnet}}, but perhaps could be useful elsewhere.
#' @examples
#' ## generate data (no intercept)
#' set.seed(4799623)
#' cn <- c()
#' for (i in 1:100) cn[i] <- paste0("x", i)
#' tb <- rbinom(100, 1, 0.05)
#' tx <- matrix(rnorm(10000), nrow = 100, ncol = 100,
#'              dimnames = list(1:100, cn))
#' ty <- tx %*% tb + rnorm(100)
#'
#' ## fit model and export model fit stats
#' cv1 <- cv.glmnet(x = tx, y = ty, family = "gaussian", alpha = 1)
#' glmnet_measure(cv1, x = tx, y = ty)
#' @export
glmnet_measure <- function(cv.model, x = NULL, y = NULL,
                           lambda.criteria = "lambda.min",
                           type.measure = c("default", "mse", "deviance",
                                            "class", "auc", "mae", "C"),
                           family = "gaussian", alpha = 1) {
  # fit the original model
  if (lambda.criteria == "lambda.min") {
    lambda = cv.model$lambda.min
  }

  if (lambda.criteria == "lambda.1se") {
    lambda = cv.model$lambda.1se
  }
  original.fit <- glmnet::glmnet(x = x, y = y, family = family,
                                 alpha = alpha, lambda = lambda)

  # model deviance
  deviance <- stats::deviance(original.fit)

  # calculate predicted values for observed design matrix
  yhat <- stats::predict(cv.model, newx = x, s = lambda.criteria)

  # mean square error
  mse <- mean((y - yhat) ^ 2)

  # mean absolute error
  mae <- mean(abs(y - yhat))

  m.df <- data.frame(model = "glmnet", s0 = lambda,
                         deviance = deviance, mse = mse, mae = mae)

  if (family == "binomial") {
    # auc
    yhatn <- as.numeric(yhat)
    yn <- as.numeric(y)
    AUC <- suppressMessages(pROC::auc(yn, yhatn))
    m.df$auc <- as.numeric(AUC)

    # misclassification
    ae <- abs(y - yhat)
    aei <- rep(0, length(y))
    aei[ae > 0.5] <- 1
    m.df$misclassification <- mean(aei)
  }

  return(m.df)
}
