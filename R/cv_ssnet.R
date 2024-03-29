#' Cross Validation for ssnet Models
#'
#' Perform k-fold cross validation for spike-and-slab elastic net models.
#'
#' @param s0,s1 A vector of user-selected possible values for the spike scale
#' and slab scale parameter, respectively. The default is
#' \code{s0 = seq(0.01, 0.1, 0.01)} and \code{s1 = 1}. However, the user
#' should select values informed by the practical context of the analysis.
#' @param nfolds Numeric value indicating the number of folds to create.
#' @param ncv Numeric value indicating the number of times to perform cross
#' validation.
#' @param foldid An (optional) vector of values between 1 and \code{nfold}
#' identifying the fold for each observation. When supplied \code{nfolds} may
#' be omitted. If \code{ncv > 1}, then supply a matrix or data frame  where
#' each column contains fold identifiers. If \code{foldid} is supplied, it
#' supersedes \code{ncv} and \code{nfolds}.
#' @param fold.seed A scalar seed value for cross validation; ensures the
#' folds are the same upon re-running the function. Alternatively, use
#' \code{foldid} to manually specify folds.
#' @inheritParams validate_ssnet
#' @inheritParams ssnet
#' @return Either a data frame of model fitness measures or a list whose
#' elements are data frames of model fitness measures and parameter estimates,
#' respectively, depending on the value of output_param_ets.
#' @examples
#' xtr <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' xte <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' b <- rnorm(5)
#'
#' ## continuous outcome
#' ytr <- xtr %*% b + rnorm(100)
#' yte <- xte %*% b + rnorm(100)
#'
#' ## binary outcome
#' ybtr <- ifelse(ytr > 0, 1, 0)
#' ybte <- ifelse(yte > 0, 1, 0)
#'
#' ## multinomial outcome
#' ymtr <- dplyr::case_when(
#'   ytr > 1 ~ "a",
#'   ytr <= 1 & ytr > -1 ~ "b",
#'   ytr <= -1 ~ "c"
#' )
#' ymte <- dplyr::case_when(
#'   yte > 1 ~ "a",
#'   yte <= 1 & yte > -1 ~ "b",
#'   yte <= -1 ~ "c"
#' )
#'
#' cv_ssnet(
#'   model = "ss", family = "gaussian",
#'   x = rbind(xtr, xte), y = c(ytr, yte),
#'   s0 = c(0.01, 0.05, 0.10), s1 = c(1, 2.5),
#'   nfolds = 3, ncv = 2
#' )
#'
#' \dontrun{
#' cv_ssnet(
#'   model = "ss", family = "binomial",
#'   x = rbind(xtr, xte), y = c(ybtr, ybte),
#'   s0 = c(0.01, 0.05, 0.10), s1 = c(1, 2.5),
#'   nfolds = 3, ncv = 2, classify = TRUE,
#'   output_param_est = TRUE
#' )
#'
#' cv_ssnet(
#'   model = "ss", family = "multinomial",
#'   x = rbind(xtr, xte), y = c(ymtr, ymte),
#'   s0 = c(0.01, 0.05, 0.10), s1 = c(1, 2.5),
#'   nfolds = 3, ncv = 2, classify = FALSE,
#'   output_param_est = TRUE
#' )
#'}
#'
#' @export
cv_ssnet <- function(
  model,
  alpha = c(0.5, 1),
  s0 = seq(0.01, 0.1, 0.01),
  s1 = c(1, 2.5),
  classify = FALSE,
  classify.rule = 0.5,
  nfolds = 10, ncv = 1,
  foldid = NULL,
  fold.seed = NULL,
  x, y, family,
  offset = NULL,
  epsilon = 1e-04,
  maxit = 50,
  init = NULL,
  group = NULL,
  Warning = FALSE,
  verbose = FALSE,
  opt.algorithm = "LBFGS",
  iar.data = NULL,
  iar.prior = FALSE,
  p.bound = c(0.01, 0.99),
  tau.prior = "none",
  stan_manual = NULL,
  lambda.criteria = "lambda.min",
  output_param_est = FALSE,
  type.multinomial = "grouped"
){
  # reproducibility
  if (is.null(fold.seed) == FALSE) {
    set.seed(fold.seed)
  }
  # for storing model fitness output
  model_fit <- NULL

  # for storing parameter estimates
  if (output_param_est == TRUE) {
    param_est <- NULL
  }

  # ensure foldid is valid or build foldid from scratch
  if (is.null(foldid) == FALSE) {
    if (
      is.numeric(foldid) == FALSE &
      is.matrix(foldid) == FALSE &
      is.data.frame(foldid) == FALSE) {
        stop("foldid must be a numeric vector, matrix, or data frame.")
        #stop("foldid must be NULL or a (numeric) vector, matrix, or data frame.")
    }
    if (is.vector(foldid) == TRUE) {
      ncv <- 1
      nfolds <- length(unique(foldid))
      foldid <- matrix(foldid, ncol = 1)
    } else {
      ncv <- ncol(foldid)
      nfolds.b <- c()
      for (nf in seq_len(ncv)) {
        nfolds.b[nf] <- length(unique(foldid[ , nf]))
      }
      if (any(nfolds.b != nfolds.b[1])) {
        stop("Each column of foldid must have the same number of unique folds.")
      } else {
        nfolds <- nfolds.b[1]
      }
    }
  } else {
    foldid <- matrix(
      sample(seq_len(nfolds),
             size = nrow(x) * ncv,
             replace = TRUE),
      ncol = ncv,
      nrow = nrow(x)
    )
  }
  # begin cv procedure
  for (i in seq_len(ncv)) {
    for (j in seq_len(nfolds)) {
      # divide data into:
      # training set (all but fold j)
      xtr.ij <- x[foldid[ , i] != j, ]
      ytr.ij <- y[foldid[ , i] != j]
      # test set (fold j)
      xte.ij <- x[foldid[ , i] == j, ]
      yte.ij <- y[foldid[ , i] == j]

      for (a in seq_len(length(alpha))) {
        for (s0b in seq_len(length(s0))) {
          if (model == "glmnet") {
            mod.ijas <- validate_ssnet(
              model = model, family = family,
              alpha = alpha[a], s0 = s0[s0b], s1 = s0[s0b],
              x.train = xtr.ij, y.train = ytr.ij,
              x.test = xte.ij, y.test = yte.ij,
              classify = classify, classify.rule = classify.rule,
              type.multinomial = type.multinomial,
              offset = offset, epsilon = epsilon,
              maxit = maxit, init = init, group = group,
              Warning = Warning, verbose = verbose,
              opt.algorithm = opt.algorithm,
              iar.data = iar.data, iar.prior = iar.prior,
              p.bound = p.bound, tau.prior = tau.prior,
              stan_manual = stan_manual, lambda.criteria = lambda.criteria,
              output_param_est = output_param_est
            )
            if (output_param_est == FALSE) {
              mod.fit.ijas <- cbind(
                model = model,
                ncv.id = i,
                fold.id = j,
                alpha = alpha[a],
                s0 = s0[s0b],
                s1 = s0[s0b],
                mod.ijas
              )
            } else {
              mod.fit.ijas <- cbind(
                model = model,
                ncv.id = i,
                fold.id = j,
                alpha = alpha[a],
                s0 = s0[s0b],
                s1 = s0[s0b],
                mod.ijas$model_fitness
              )
              param.est.ijas <- cbind(
                model = model,
                ncv.id = i,
                fold.id = j,
                alpha = alpha[a],
                s0 = s0[s0b],
                s1 = s0[s0b],
                mod.ijas$param_est
              )

              param_est <- dplyr::bind_rows(
                param_est,
                param.est.ijas
              )
            }
            model_fit <- dplyr::bind_rows(
              model_fit,
              mod.fit.ijas
            )
          } else {
            for (s1b in seq_len(length(s1))) {
              mod.ijas <- validate_ssnet(
                model = model, family = family,
                alpha = alpha[a], s0 = s0[s0b], s1 = s1[s1b],
                x.train = xtr.ij, y.train = ytr.ij,
                x.test = xte.ij, y.test = yte.ij,
                classify = classify, classify.rule = classify.rule,
                type.multinomial = type.multinomial,
                offset = offset, epsilon = epsilon,
                maxit = maxit, init = init, group = group,
                Warning = Warning, verbose = verbose,
                opt.algorithm = opt.algorithm,
                iar.data = iar.data, iar.prior = iar.prior,
                p.bound = p.bound, tau.prior = tau.prior,
                stan_manual = stan_manual,
                lambda.criteria = lambda.criteria,
                output_param_est = output_param_est
              )
              if (output_param_est == FALSE) {
                mod.fit.ijas <- cbind(
                  model = model,
                  ncv.id = i,
                  fold.id = j,
                  alpha = alpha[a],
                  s0 = s0[s0b],
                  s1 = s1[s1b],
                  mod.ijas
                )
              } else {
                mod.fit.ijas <- cbind(
                  model = model,
                  ncv.id = i,
                  fold.id = j,
                  alpha = alpha[a],
                  s0 = s0[s0b],
                  s1 = s1[s1b],
                  mod.ijas$model_fitness
                )
                param.est.ijas <- cbind(
                  model = model,
                  ncv.id = i,
                  fold.id = j,
                  alpha = alpha[a],
                  s0 = s0[s0b],
                  s1 = s1[s1b],
                  mod.ijas$param_est
                )

                param_est <- dplyr::bind_rows(
                  param_est,
                  param.est.ijas
                )
              }
              model_fit <- dplyr::bind_rows(
                model_fit,
                mod.fit.ijas
              )
            }
          }
        }
      }
    }
  }
  if (output_param_est == TRUE) {
    return(
      list(
        model_fit = model_fit,
        param_est = param_est
      )
    )
  } else {
    return(model_fit)
  }
}
