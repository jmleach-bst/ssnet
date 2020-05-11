#' Select the Optimal Model from Each Simulation
#'
#' For a set of simulations that used multiple spike scale parameters this function selects the optimal
#' model based on user-selected criteria and creates a new data frame of statistics for the optimal models.
#' If desired, outputs summaries from the new data frame.
#'
#' @param sim.data A data frame containing simulation results. The data is intended to be generated
#'  from stacking multiple \code{compare_ssnet()} outputs in a single data frame.
#' @param criteria Specifies the criteria for model selection. Options are \code{"deviance"}, \code{"mse"},
#' \code{"mae"} for deviance, mean-square error, and mean absolute error, respectively. When
#' \code{family = "binomial"}, additional options are \code{"auc"} and \code{"misclassification"}, for
#' Area under the ROC curve and the percentage of cases where the difference between the observed and
#' predicted values is greater than 1/2.
#' @param tie.breaker This argument decides how to break ties when multiple values of spike scale
#' minimize \code{criteria}. The default is \code{"min"}, which selects the smallest spike scale that
#' minimizes \code{criteria}; alternatively, \code{"max"}, which selects the largest spike scale that
#' minimizes \code{criteria}.
#' @param print.details Logical. Determines whether to print periodic output.
#' @param suppress.warnings Logical. When \code{TRUE} suppresses warning output. Default is \code{FALSE}.
#' @return A data frame containing where each row is an optimal model from a set of models.
#' @examples
#' sdt <- data.frame(model = rep("mystery", 12),
#'                   s0 = rep(c(0.01, 0.02, 0.03), 4),
#'                   deviance = c(-0.17267380, -0.00398748, 0.49164019, -0.19134208, 0.05202775, -0.76916548, 0.08757465,
#'                                -0.73574828, -0.64348907, -0.46910645, -0.56955834, -1.54488900))
#' optimal_s0(sdt, criteria = "deviance")
#'
#' @export
optimal_s0 <- function(sim.data, criteria,
                       tie.breaker = "min",
                       print.details = FALSE,
                       suppress.warnings = FALSE) {

  # number of spike scale values fit
  n.s0 <- length(unique(sim.data$s0))

  # number of simulated datasets
  n.sim <- nrow(sim.data) / n.s0

  if (print.details == TRUE) {
    cat(n.s0, "spike scale values. \n")
    cat(n.sim, "simulated datasets. \n")
  }

  # generate simulation id's for easy looping
  sim.id <- c()
  for (i in 1:n.sim) {
    sim.id <- c(sim.id, rep(i, n.s0))
  }
  sim.data$sim.id <- sim.id

  # pick optimal model in each dataset
  for (m in 1:n.sim) {
    sim.data.m <- sim.data[sim.data$sim.id == m, ]
    mod.m <- sim.data.m[sim.data.m[criteria] == min(sim.data.m[criteria]), ]

    if (nrow(mod.m) > 1) {
      if (suppress.warnings == FALSE) {
        warning(cat("Simulation run ", m, " has multiple optimal models; selecting model with ", tie.breaker, "s0. \n"))
      }
      if (tie.breaker == "min") {
        mod.m <- mod.m[1, ]
      }
      if (tie.breaker == "max") {
        mod.m <- mod.m[nrow(mod.m), ]
      }
    }

    # if (print.details == TRUE) {
    #   cat("Optimal model for dataset ", m, " has s0 = ", mod.m$s0, "\n")
    # }
    if (print.details == TRUE & (mod.m$s0 == max(sim.data.m$s0) | mod.m$s0 == min(sim.data.m$s0))) {
      if (suppress.warnings == FALSE) {
        warning(cat("Simulation run ", m, "selected minimum or maximum allowable s0 = ", mod.m$s0, ". \n "))
        print(sim.data.m)
      }
    }
    if (m > 1) {
      out <- rbind(out, mod.m)
    } else {
      out <- mod.m
    }
  }
  return(out)
}

