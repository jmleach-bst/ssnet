#' Generate Folds for K-Fold Cross Validation
#'
#' Unchanged from \code{BhGLM}, but was not in NAMESPACE, so reproduced here.
#' @param nobs The number of observations.
#' @param foldid An array of fold identifiers. If not supplied by the user, generated within function.
#' @param fold.seed An integer that sets the seed for generating folds.
#' @inheritParams BhGLM::glmNet
#'
generate.foldid <- function (nobs, nfolds = 10, foldid = NULL, ncv = 1, fold.seed = NULL)
{
  if (nfolds > nobs)
    nfolds <- nobs
  if (nfolds == nobs)
    ncv <- 1
  if (is.null(foldid)) {
    if (is.null(fold.seed) == FALSE) {
      set.seed(fold.seed)
    }
    foldid <- array(NA, c(nobs, ncv))
    for (j in 1:ncv) foldid[, j] <- sample(rep(seq(nfolds),
                                               length = nobs))
  }
  foldid <- as.matrix(foldid)
  nfolds <- max(foldid)
  ncv <- ncol(foldid)
  list(foldid = foldid, nfolds = nfolds, ncv = ncv)
}
