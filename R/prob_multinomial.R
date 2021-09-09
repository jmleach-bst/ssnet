#' Calculate Probabilities from Multinomial Regression
#'
#' Uses design matrix and parameter estimates to calculate probabilities for each class
#' for each subject.
#'
#' @param x Design, or input, matrix, of dimension nobs x nvars; each row is an observation vector. It
#' is recommended that \code{x} have user-defined column names for ease of identifying variables. If
#' missing, then \code{colnames} are internally assigned \code{x1}, \code{x2}, ... and so forth.
#' @param b A list whose elements are vectors of parameter estimates and whose names correspond to the
#' values of \code{y}. If \code{fitted.model} is obtained from \code{ssnet}, then these estimates
#' are stored in \code{fitted.model$coefficients}.
#' @param a0 A matrix of with \code{nrow} equal the number of class and a single column containing
#' the intercept estimates for each class. If \code{fitted.model} is obtained from \code{ssnet},
#' then these estimates are stored in \code{fitted.model$a0}. Note that the row names should
#' correspond to the class labels from \code{y}.
#' @param classify Logical. When \code{TRUE} (default), appends a column containing the class with the
#' largest predicted probability.
#'
#' @return A data frame where each column contains predicted class probabilities for each subject.
#' @details The parameterization employed here is as in Friedman (2010), i.e., over-parameterized,
#' which under the elastic net penalty allows for unique parameter estimates. Thus, there is not a
#' reference class, but each class has an associated set of parameter estimates. Probabilities are
#' then obtained by \eqn{Pr(Y_i = k | X_i) = \exp(X_i^T B_k) / \sum_{v=1}^V \exp(X_i^T B_v)} where
#' X_i is a the design matrix for subject i, V is the number of classes, and B_v is the vector of
#' parameters associated with class v.
#'
#' @examples
#' x <- matrix(rnorm(10*5), nrow = 10, ncol = 5)
#' colnames(x) <- paste0("x", 1:ncol(x))
#' x0 <- cbind(x0 = 1, x)
#' b <- list(a = runif(5), b = runif(5), c = runif(5))
#' for (i in 1:length(b)) {
#'     names(b[[i]]) <- colnames(x)
#' }
#' a0 <- matrix(runif(3), nrow = 3, ncol= 1)
#' prob_multinomial(x = x, b = b, a0 = a0)
#'
#' @references
#' \insertRef{Friedman:2010}{ssnet}
#'
#' @export
prob_multinomial <- function(x, b, a0,
                             classify = TRUE) {

  # check matrix
  if (is.matrix(x) == FALSE) {
    stop("x should be a matrix")
  }

  # assign variable names if necessary
  if (is.null(colnames(x)) == TRUE) {
    colnames(x) <- paste0("x", 1:ncol(x))
  }

  # add column of 1's for intercept
  x0 <- cbind(x0 = 1, x)

  # require class names
  if (is.null(names(b))) {
    stop("b should have names corresponding to class levels.")
  }

  if (any(names(b) != rownames(a0))) {
    stop("b and a0 should have the same class level names.")
  }

  # combine intercept + parameter estimates into 1 df
  b.mat <- as.matrix(
    cbind(a0, data.frame(dplyr::bind_rows(b)))
    )
  colnames(b.mat) <- colnames(x0)
  rownames(b.mat) <- names(b)
  # print(b.mat)

  # compute & store raw predictions
  exp.xb <- list()
  for (i in 1:nrow(b.mat)) {
    exp.xb[[i]] <- exp(x0 %*% b.mat[i, ])
  }
  names(exp.xb) <- rownames(b.mat)
  # print(exp.xb)
  exp.xb.df <- data.frame(dplyr::bind_cols(exp.xb))
  # print(exp.xb.df)
  denominator <- rowSums(exp.xb.df)

  pr.yi <- list()
  for (i in 1:ncol(exp.xb.df)) {
    pr.yi[[i]] <- exp.xb.df[, i] / denominator
  }
  names(pr.yi) <- colnames(exp.xb.df)
  pr.yi.df <- data.frame(dplyr::bind_cols(pr.yi))

  if (classify == TRUE) {
    predicted.class <- c()
    for (i in 1:nrow(pr.yi.df)) {
      predicted.class[i] <- colnames(pr.yi.df)[which(
        pr.yi.df[i, ] == max(pr.yi.df[i, ])
      )]
    }
    df <- cbind(pr.yi.df, predicted.class)
    colnames(df) <- c(names(b),
                      "predicted.class")
    return(df)
  } else {
    return(pr.yi.df)
  }

}
























