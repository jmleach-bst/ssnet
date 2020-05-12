#' Plot Probabilities A Model Keeps Parameters
#'
#' Provides a plot to visualize the probability that each of several models includes specific parameters.
#'
#' @param data A data frame whose first column should be named \code{"model"} and contain model names, and whose
#' remaining columns contain the probabilities that the models included that parameter.
#' @param grp.col A vector of color names that will correspond to each model.
#' @param remove.int Logical. When \code{TRUE} assumes that the second column is for the intercept and removes it.
#' The default is \code{FALSE}.
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_x_continuous scale_y_continuous ggtitle theme
#' @return An object of class gg and ggplot.
#' @examples
#' dfp <- data.frame(model = c("glmnet", "bmlasso", "bmlasso_iar"),
#'                   x1 = c(0.01, 0.05, 0.10),
#'                   x2 = c(0.05, 0.10, 0.01),
#'                   x3 = c(0.95, 0.90, 0.85),
#'                   x4 = c(0.7, 0.6, 0.5),
#'                   x5 = c(0.1, 0.2, 0.3))
#' dfp$model <- factor(dfp$model, levels = c("glmnet", "bmlasso", "bmlasso_iar"))
#' plot_param(data = dfp, grp.col = c("#E69F00", "#56B4E9", "#009E73"))
#'
#'@export
plot_param <- function(data, grp.col, remove.int = FALSE) {

  # remove intercept
  if (remove.int == TRUE) {
    data[, -2]
  }

  # rearrage the data
  mdfp <- as.matrix(data[, -1])
  prob <- as.vector(mdfp)
  mi <- expand.grid(unique(data$model), 1:ncol(mdfp))
  dfp2 <- data.frame(cbind(prob, mi))
  names(dfp2) <- c("Probability", "Model", "x.index")

  # this is silly, but otherwise "no visible binding" warning
  Probability <- dfp2$Probability
  x.index <- dfp2$x.index
  Model <- dfp2$Model

  # obtain the plots
  ggplot2::ggplot(data = dfp2,
                  mapping = ggplot2::aes(y = Probability,
                                         x = x.index,
                                         color = Model)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = grp.col) +
    ggplot2::scale_x_continuous(name = "Parameter Index") +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    ggplot2::ggtitle("Probability the Parameter is Left in the Model") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), text = ggplot2::element_text(size = 12))
}
