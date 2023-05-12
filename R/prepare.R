#' Prepare the model for analysis.
#'
#' @param x The design, or input, matrix, of dimension nobs x nvars; each row
#' is an observation vector.
#' @param intercept Logical. When \code{intercept = TRUE}, an intercept is
#' included in the model.
#' @param prior.mean A vector of prior means for the parameters. If
#' \code{length(prior.mean) = 1}, then all parameters are assumed to have a
#' common prior mean.
#' @param prior.sd A vector defining prior standard deviations in the normal
#' priors of the coefficients. If provided, they are starting values of the
#' prior standard deviations for the iterative algorithms. This argument is
#' used within the function \code{bglm()}.If \code{length(prior.sd) = 1}, then
#' all parameters are assumed to have a common prior standard deviation.
#' @param prior.scale A vector containing the prior scale values for double
#' exponential or t prior. In the spike-and-slab lasso model, this is
#' generally only the slab prior scale.
#' @param prior.df A scalar defining the prior degrees of freedom when the
#' t-distribution is used as a prior for the parameters.
#' @param multinomial Logical. When \code{TRUE}, prepare the data for analysis
#' of multinomial outcome according to parameterization in \code{glmnet}, i.e.,
#' a unique set of parameters is associated with each outcome category. See
#' section 4 of \insertCite{Friedman:2010}{ssnet} for details.
#' @param outcome.cats Either an integer denoting the number of response
#' categories, or a vector containing outcome categories. Only used when
#' \code{multinomial = TRUE}.
#' @param group A numeric vector, or an integer, or a list indicating the
#' groups of predictors. If \code{group = NULL}, all the predictors form a
#' single group. If \code{group = K}, the predictors are evenly divided into
#' groups each with K predictors. If group is a numberic vector, it defines
#' groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2:
#' \code{(group[2]+1):group[3]}, Group 3: \code{(group[3]+1):group[4]}, ...
#' If group is a list of variable names, \code{group[[k]]} includes variables
#' in the k-th group. The mixture double-exponential prior is only used for
#' grouped predictors. For ungrouped predictors, the prior is
#' double-exponential with scale \code{ss[2]} and mean 0.
#' @return A list containing necessary elements for Bayesian analysis in the
#' package \code{BhGLM}.
#' @note This function is a modified version the function \code{prepare} from
#' the R package \code{BhGLM}.
#' @references
#'
#' \insertRef{Friedman:2010}{ssnet}
#'
#' @export
prepare <- function(
    x,
    intercept,
    prior.mean,
    prior.sd,
    prior.scale,
    prior.df,
    group,
    multinomial = FALSE,
    outcome.cats = NULL
){
  x0 <- x
  if (intercept == TRUE) {
    x0 <- x[, -1, drop = FALSE]
  }

  if (is.null(colnames(x0)) == TRUE) {
    internal.xnames <- c()
    for (i in seq_len(ncol(x0))) {
      internal.xnames[i] <- paste0("x", i)
    }
    colnames(x0) <- internal.xnames
  }

  # handle grouping tasks
  g <- Grouping(all.var = colnames(x0), group = group)
    group <- g$group
    group.vars <- g$group.vars
    ungroup.vars <- g$ungroup.vars
    covars <- g$ungroup.vars
    if (is.list(group)) {
      if (length(unlist(group)) > length(unique(unlist(group)))) {
        x1 <- as.data.frame(x0)
        x1 <- x1[, c(covars, unlist(group))]
        g <- c(length(ungroup.vars), length(ungroup.vars) +
                 cumsum(lapply(group, length)))
        for (j in seq_len(length(group) - 1)) {
          group.vars[[j]] <- colnames(x1[, (g[j] + 1):g[j + 1]])
        }
        x1 <- as.matrix(x1)
        x <- x1
        if (intercept == TRUE) {
          x <- cbind(1, x)
          colnames(x)[1] <- "(Intercept)"
        }
      }
    }

  # priors for model parameters
  J <- ncol(x)
  if (intercept == TRUE & J > 1) {
    prior.mean <- c(0, prior.mean)
    prior.scale <- c(prior.scale[1], prior.scale)
    prior.df <- c(prior.df[1], prior.df)
  }
  if (length(prior.mean) < J) {
    prior.mean <- c(
      prior.mean,
      rep(prior.mean[length(prior.mean)],
          J - length(prior.mean)
          )
      )
  }

  if (length(prior.scale) < J) {
    prior.scale <- c(
      prior.scale,
      rep(prior.scale[length(prior.scale)],
          J - length(prior.scale)
          )
      )
  }

  if (length(prior.df) < J) {
    prior.df <- c(
      prior.df,
      rep(prior.df[length(prior.df)],
          J - length(prior.df)
          )
      )
  }

  prior.mean <- prior.mean[seq_len(J)]
  prior.scale <- prior.scale[seq_len(J)]
  prior.df <- prior.df[seq_len(J)]
  prior.df <- ifelse(prior.df == Inf, 1e+10, prior.df)

  if (is.null(prior.sd)) {
    prior.sd <- prior.scale + 0.2
  }

  if (length(prior.sd) < J) {
    prior.sd <- c(
      prior.sd,
      rep(prior.sd[length(prior.sd)],
          J - length(prior.sd)
          )
      )
  }

  prior.sd <- prior.sd[seq_len(J)]
  sd.x <- apply(x, 2, sd, na.rm = TRUE)
  min.x.sd <- 1e-04
  prior.sd <- ifelse(sd.x < min.x.sd, 1e-04, prior.sd)

  if (intercept == TRUE) {
    prior.sd[1] <- 1e+10
  }

  names(prior.mean) <- names(prior.scale) <- names(prior.df) <- names(prior.sd) <- colnames(x)
  if (intercept == TRUE) {
    covars <- c(colnames(x)[1], covars)
  }

  if (!is.null(covars)) {
    prior.mean[covars] <- 0
  }

  if (multinomial == TRUE) {
    if (is.integer(outcome.cats) & length(outcome.cats) == 1) {
      num.cats <- outcome.cats
    } else {
      num.cats <- length(unique(outcome.cats))
    }

    prior.mean.list <- list()
    prior.sd.list <- list()
    for (i in seq_len(num.cats)) {
      prior.mean.list[[i]] <- prior.mean
      prior.sd.list[[i]] <- prior.sd
    }

    return(
      list(
        x = x,
        prior.mean = prior.mean.list,
        prior.sd = prior.sd.list,
        prior.scale = prior.scale,
        prior.df = prior.df,
        sd.x = sd.x,
        min.x.sd = min.x.sd,
        group = group,
        group.vars = group.vars,
        ungroup.vars = ungroup.vars
        )
      )
  } else {
    return(
      list(
        x = x,
        prior.mean = prior.mean,
        prior.sd = prior.sd,
        prior.scale = prior.scale,
        prior.df = prior.df,
        sd.x = sd.x,
        min.x.sd = min.x.sd,
        group = group,
        group.vars = group.vars,
        ungroup.vars = ungroup.vars
        )
      )
  }
}
