#' Prepare variables for group-based inference.
#'
#' This function uses the user-defined groups to properly prepare the variables for
#' group-specific estimation/inference.
#'
#' @param all.var A vector containing the names of the variables in the design matrix.
#' @param group A numeric vector, or an integer, or a list indicating the groups of predictors.
#' If \code{group = NULL}, all the predictors form a single group. If \code{group = K}, the
#' predictors are evenly divided into groups each with K predictors. If group is a numberic vector,
#' it defines groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2: \code{(group[2]+1):group[3]},
#' Group 3: \code{(group[3]+1):group[4]}, ... If group is a list of variable names, \code{group[[k]]}
#' includes variables in the k-th group. The mixture double-exponential prior is only used for grouped
#' predictors. For ungrouped predictors, the prior is double-exponential with scale \code{ss[2]} and mean 0.
#' @return A list.
#' @note This function is taken from the R package \code{BhGLM}.
Grouping <- function (all.var, group)
{
  n.vars <- length(all.var)
  group.vars <- list()
  if (is.list(group))
    group.vars <- group
  else {
    if (is.numeric(group) & length(group) > 1) {
      group <- sort(group)
      if (group[length(group)] > n.vars)
        stop("wrong grouping")
    }
    if (is.numeric(group) & length(group) == 1)
      group <- as.integer(seq(0, n.vars, length.out = n.vars / group + 1))
    if (is.null(group))
      group <- c(0, n.vars)
    group <- unique(group)
    for (j in 1:(length(group) - 1)) group.vars[[j]] <- all.var[(group[j] +
                                                                   1):group[j + 1]]
  }
  all.group.vars <- unique(unlist(group.vars))
  if (length(all.group.vars) == n.vars)
    ungroup.vars <- NULL
  else ungroup.vars <- all.var[which(!all.var %in% all.group.vars)]
  group.new <- c(length(ungroup.vars), length(ungroup.vars) +
                   cumsum(lapply(group.vars, length)))
  var.new <- c(ungroup.vars, unlist(group.vars))
  list(group = group, group.vars = group.vars, ungroup.vars = ungroup.vars,
       group.new = group.new, var.new = var.new)
}
