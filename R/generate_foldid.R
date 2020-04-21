#' Generate Folds for K-Fold Cross Validation
#'
#' Unchanged from \code{BhGLM}, but was not in NAMESPACE, so reproduced here.
#'
generate.foldid <- function (nobs, nfolds = 10, foldid = NULL, ncv = 1)
{
  if (nfolds > nobs)
    nfolds <- n
  if (nfolds == nobs)
    ncv <- 1
  if (is.null(foldid)) {
    foldid <- array(NA, c(nobs, ncv))
    for (j in 1:ncv) foldid[, j] <- sample(rep(seq(nfolds),
                                               length = nobs))
  }
  foldid <- as.matrix(foldid)
  nfolds <- max(foldid)
  ncv <- ncol(foldid)
  list(foldid = foldid, nfolds = nfolds, ncv = ncv)
}
