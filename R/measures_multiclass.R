#' Obtain Measures of Model Performance for Mutliclass Classification
#'
#' Compare observed and predicted classes in cases where there are more than 2 classes.
#'
#' @param y Scalar response variable denoting class membership. It is preferred that \code{y} is a
#' factor with two or more levels, but will modify internally for numeric or character variables.
#' @param pr_yi A data frame or matrix where each row is a subject observation where columns indicate
#' the probability the subjects belongs to the respective class. The column names should match the
#' factor levels for \code{y}. This information can be obtained using \code{prob_multinomial()}.
#' @param y_hat A vector of predicted classes that should have the same factor levels as \code{y}.
#' When \code{y_hat} is supplied it supercedes \code{pr_yi} for classification. However, \code{pr_yi}
#' will still be used to calculate the deviance.
#' @param print_check Logical. When \code{TRUE}, prints intermediate results.
#' @note While this function will work when there are only 2 classes, it is recommended to use
#' \code{eval_classify()} when there are only 2 classes.
#' @return A data frame with a single row containing columns for average accuracy (avg_acc), average
#' per-class classification error (pce), micro-averaged positive predictive value (ppv_micro),
#' sensitivity (sn_micro), F1 score (f1_micro), and macro-averaged positive predictive value (ppv_macro),
#' sensitivity (sn_macro), F1 score (f1_macro).
#' @details Multiclass measures are calculated according to Sokolova & Lapalme (2009). We give brief
#' details here, but refer users to Table 3 of Sokolova & Lapalme (2009) for precise definitions.
#' @examples
#'
#' n <- 500
#' y <- sample(c("ack", "eek", "ahh"), size = n, replace = TRUE)
#' x <- matrix(rnorm(n*5), nrow = n, ncol = 5)
#' colnames(x) <- paste0("x", 1:ncol(x))
#' mfit <- glmnet::glmnet(x = x, y = y, family = "multinomial",
#'                        type.multinomial = "grouped",
#'                        lambda = 0.01)
#' xnew <- matrix(rnorm(n*5), nrow = n, ncol = 5)
#' b <- list()
#' for (i in 1:length(mfit$beta)) {
#'   b[[i]] <- as.numeric(mfit$beta[[i]])
#' }
#' names(b) <- names(mfit$beta)
#' for (i in 1:length(b)) {
#'     names(b[[i]]) <- colnames(x)
#' }
#' pm <- prob_multinomial(x = xnew, b = b, a0 = mfit$a0)
#' measures_multiclass(y = y, pr_yi = pm[,-4])
#' measures_multiclass(y = y, y_hat = pm$predicted.class)
#'
#' @references
#' \insertRef{Sokolova+Lapalme:2009}{ssnet}
#' @export
measures_multiclass <- function(y, pr_yi = NULL, y_hat = NULL,
                                print_check = FALSE) {

  ##########
  # Checks #
  ##########

  # convert y to factor
  if (is.factor(y) == FALSE) {
    y <- as.factor(y)
  }

  if (length(levels(y)) == 2){
    warning("Recommended to use eval_classify() when only 2 classes.")
  }

  if (is.null(pr_yi) == TRUE & is.null(y_hat) == TRUE) {
    stop("User must supply either pr_yi or y_hat.")
  }

  if (is.null(y_hat) == FALSE) {
    # y_hat supercedes pr_yi
    if (is.null(pr_yi) == FALSE) {
      warning("Recommended to supply either y_hat or pr_yi. Using y_hat and ignoring pr_yi.")
      pr_yi <- NULL
    }

    if (any(levels(y) != levels(y_hat))) {
      stop("Factor levels of y_hat should correspond to factor levels of y.")
    }

  }

  if (is.null(pr_yi) == FALSE) {
    if (any(levels(y) != colnames(pr_yi))) {
      stop("Column names for predicted probabilties should correspond to factor levels of y.")
    }

    y_hat <- c()
    for (i in 1:nrow(pr_yi)) {
      pc.i <- which(pr_yi[i, ] == max(pr_yi[i, ]))
      if (length(pc.i) > 1) {
        print(pr_yi[i, ])
        pc.i <- sample(pc.i, size = 1, replace = FALSE)
        warning("Multiple categories have maximum probability. Class selected randomly.")
      }
      y_hat[i] <- colnames(pr_yi)[pc.i]
    }
  }

  if (length(y) != length(y_hat)) {
    stop("Observed y and predicted y_hat should be vectors of equal length.")
  }

  # convert y to factor
  if (is.factor(y_hat) == FALSE) {
    y_hat <- as.factor(y_hat)
  }

  if (print_check == TRUE) {
    print(table(y_hat))
  }

  #############################################
  # Calculating class-specific TP, FP, TN, FN #
  #############################################

  # confusion.matrices <- list()
  out_pn <- list()
  for (i in 1:length(levels(y))) {
    y_i <- ifelse(y == levels(y)[i], 1, 0)
    y_hat_i <- ifelse(y_hat == levels(y)[i], 1, 0)
    tp_i <- sum(ifelse(y_i == 1 & y_hat_i == 1, 1, 0))
    fp_i <- sum(ifelse(y_i == 0 & y_hat_i == 1, 1, 0))
    tn_i <- sum(ifelse(y_i == 0 & y_hat_i == 0, 1, 0))
    fn_i <- sum(ifelse(y_i == 1 & y_hat_i == 0, 1, 0))

    if (print_check == TRUE) {
      print(cbind(y, y_hat, y_i, y_hat_i,
                  tp = ifelse(y_i == 1 & y_hat_i == 1, 1, 0),
                  fp = ifelse(y_i == 0 & y_hat_i == 1, 1, 0),
                  tn = ifelse(y_i == 0 & y_hat_i == 0, 1, 0),
                  fn = ifelse(y_i == 1 & y_hat_i == 0, 1, 0)))
    }

    # cm_i <- matrix(c(tp_i, fn_i,
    #                  fp_i, tn_i),
    #                byrow = TRUE,
    #                nrow = 2,
    #                ncol = 2)
    out_pn[[i]] <- data.frame(
      tp = tp_i,
      fp = fp_i,
      tn = tn_i,
      fn = fn_i)
  }

  names(out_pn) <- levels(y)
  out_pn_df <- data.frame(dplyr::bind_rows(out_pn))
  rownames(out_pn_df) <- levels(y)

  #########################
  # Measurement Functions #
  #########################

  accuracy <- function(tp, tn, fp, fn) {
    (tp + tn) / (tp + fn + fp + tn)
  }

  error_rate <- function(tp, tn, fp, fn) {
    (fp + fn) / (tp + fn + fp + tn)
  }

  ppv <- function(tp, fp) {
    tp / (tp + fp)
  }

  sensitivity <- function(tp, fn) {
    tp / (tp + fn)
  }

  f1 <- function(tp, fn, fp) {
    (2 * tp) / (2 * tp + fn + fp)
  }

  f1_v2 <- function(ppv, sn) {
    (2 * ppv * sn) / (ppv + sn)
  }

  measurements <- list()
  for (i in 1:length(levels(y))) {
    pn_i <- out_pn_df[i, ]
    measurements[[i]] <- data.frame(
      accuracy = accuracy(
        tp = pn_i$tp,
        tn = pn_i$tn,
        fp = pn_i$fp,
        fn = pn_i$fn
      ),
      error_rate = error_rate(
        tp = pn_i$tp,
        tn = pn_i$tn,
        fp = pn_i$fp,
        fn = pn_i$fn
      ),
      ppv = ppv(
        tp = pn_i$tp,
        fp = pn_i$fp
      ),
      sensitivity = sensitivity(
        tp = pn_i$tp,
        fn = pn_i$fn
      )
    )
  }
  names(measurements) <- levels(y)
  measurements_df <- data.frame(dplyr::bind_rows(measurements))
  rownames(measurements_df) <- levels(y)
  # print(measurements_df)
  macro <- data.frame(t(apply(measurements_df, 2, mean)))

  avg_acc <- macro$accuracy
  pce <- macro$error_rate
  ppv_macro <- macro$ppv
  sn_macro <- macro$sensitivity
  f1_macro <- f1_v2(ppv = ppv_macro, sn = sn_macro)

  micro <- data.frame(t(apply(out_pn_df, 2, sum)))
  ppv_micro <- micro$tp / (micro$tp + micro$fp)
  sn_micro <- micro$tp / (micro$tp + micro$fn)
  f1_micro <- f1_v2(ppv = ppv_micro, sn = sn_micro)

  if (print_check == TRUE) {
    print(out_pn_df)
    print(measurements_df)
    print(macro)
    print(micro)
  }

  # calculate multinomial deviance
  if (is.null(pr_yi) == FALSE) {
    lp <- c()
    for (i in 1:nrow(pr_yi)) {
      lp[i] <- pr_yi[i, y[i]]
    }
    logL <- log(lp)
    deviance = -2 * sum(logL)
    return(
      data.frame(
        deviance = deviance,
        avg_acc = avg_acc,
        pce = pce,
        ppv_macro = ppv_macro,
        sn_macro = sn_macro,
        f1_macro = f1_macro,
        ppv_micro = ppv_micro,
        sn_micro = sn_micro,
        f1_micro = f1_micro
      )
    )
  } else {
    return(
      data.frame(
        avg_acc = avg_acc,
        pce = pce,
        ppv_macro = ppv_macro,
        sn_macro = sn_macro,
        f1_macro = f1_macro,
        ppv_micro = ppv_micro,
        sn_micro = sn_micro,
        f1_micro = f1_micro
      )
    )
  }
}
