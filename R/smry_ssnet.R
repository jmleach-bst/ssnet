#' Extract Simulation Summaries
#'
#' Using the tidyverse to summarize simulation results.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map map_df
#' @importFrom tidyr nest
#' @importFrom dplyr group_by mutate
#' @param sim.data A data frame containing simulation resuls. The data is
#' intended to be generated from stacking multiple \code{compare_ssnet()}
#' outputs in a single data frame.
#' @param output When \code{"raw"} outputs the nested data frame along with
#' data frames contiaining means and standard deviations of model fit
#' statistics. When \code{"mean"} or \code{"sd"} outputs a data frame of the
#' means or standard deviations, respectively, of model fit statistics.
#' Defaults to \code{"raw"}.
#' @return A data frame.
#' @examples
#' ## number of simulations
#' M <- 2
#' ## subj/sim
#' N <- 30
#' ## image resolution for spatial predictors
#' ir <- c(5, 5)
#' L <- sim2Dpredictr::chol_s2Dp(im.res = ir, rho = 0.90,
#'                               corr.structure = "ar1",
#'                               triangle = "lower")
#'
#' ## generate non-zero parameters with spatial clustering
#' betas <- sim2Dpredictr::beta_builder(index.type = "ellipse",
#'                                      w = 2, h = 2,
#'                                      row.index = 3, col.index = 3,
#'                                      B.values = 0.5, im.res = ir)
#' ## generate data
#' set.seed(68741)
#' for (m in 1:M) {
#'   datm <- sim2Dpredictr::sim_Y_MVN_X(N = N, dist = "gaussian",
#'                                      L = L$L, S = L$S,
#'                                      B = betas$B)
#'
#'   mod.out.m <- compare_ssnet(
#'     x = as.matrix(datm[, grep("X.*", names(datm), perl = TRUE)]),
#'     y = datm$Y, models = c("glmnet", "ss"),
#'    s0 = seq(0.01, 0.05), nfolds = 3,
#'    family = "gaussian", model_fit = "all",
#'    variable_selection = TRUE,
#'    B = betas$B[-1],
#'    iar.data = model_info
#'    )
#'   if (m > 1) {
#'     mod.out <- rbind(mod.out, mod.out.m)
#'   } else {
#'     mod.out <- mod.out.m
#'   }
#' }
#'
#' ## summarize measures of model fitness
#' smry_ssnet(mod.out)
#'
#' ## just means
#' smry_ssnet(mod.out, output = "mean")
#' ## just sd
#' smry_ssnet(mod.out, output = "sd")
#' @export
smry_ssnet <- function(
    sim.data,
    output = "raw"
){

  # to get around the 'no visible binding' problem
  model <- NULL
  s0 <- NULL
  data <- NULL

  # glmnet is easier to summarize b/c we don't have
  # multiple scale values for each dataset.
  if ("glmnet" %in% sim.data$model) {
    m.lamb <- mean(sim.data$s0[sim.data$model == "glmnet"])
    sd.lamb <- sd(sim.data$s0[sim.data$model == "glmnet"])
    cat(
      "Note: glmnet labmda penalty had mean = ", m.lamb,
      " and sd = ", sd.lamb, "\n")
    sim.data$s0[sim.data$model == "glmnet"] <- m.lamb
  }
  # group by model and spike scale, then nest
  ssnet.sim <- sim.data %>%
      dplyr::group_by(model, s0) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        sim.mean = purrr::map(
          data,
          .f = function(x) purrr::map_df(x, .f = mean)),
        sim.sd = purrr::map(
          data,
          .f = function(x) purrr::map_df(x, .f = sd))
        )

  if (output == "raw") {
    return(ssnet.sim)
  } else {
    if (output == "mean") {
      extr.mean <- purrr::map_dfr(ssnet.sim$sim.mean, rbind)
      out.mean <- cbind(
        model = ssnet.sim$model,
        s0 = ssnet.sim$s0,
        extr.mean
        )
      return(out.mean)
    }
    if (output == "sd") {
      extr.sd <- purrr::map_dfr(ssnet.sim$sim.sd, rbind)
      out.sd <- cbind(
        model = ssnet.sim$model,
        s0 = ssnet.sim$s0,
        extr.sd
      )
      return(out.sd)
    }
  }
}

