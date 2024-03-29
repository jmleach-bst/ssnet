#' Train and evaluate models
#'
#' Fit a model using ssnet with training data and evaluate using test data.
#'
#' @importFrom glmnet glmnet
#' @param model Specify which model to fit. Options include
#' \code{c("glmnet", "ss", "ss_iar")}.
#' @param s0,s1 A numeric value. When fitting a spike-and-slab model,
#' \code{s0} is the spike scale and \code{s1} is the slab scale. Default is
#' \code{s0 = 0.01} and \code{s1 = 1}. When \code{model = "glmnet"}, only
#'  \code{s0} is used.
#' @param x.train,x.test Design matrices for training and test data,
#' respectively.
#' @param y.train,y.test Response/outcome vectors for training and testing,
#' respectively.
#' @param output_param_est Logical. When \code{TRUE} adds an element to the
#' output list that includes parameter estimates for the fitted model.
#' Defaults is \code{FALSE}.
#' @param output_probs Logical. When \code{TRUE} and
#' \code{family = "multinomial"} adds an element to the output list that
#' contains the probabilties of being a member of each category for each
#' subject, in addition to their classification. Default is \code{FALSE}.
#' @inheritParams ssnet
#' @inheritParams glmnet::glmnet
#' @inheritParams glmnet::cv.glmnet
#' @inheritParams BhGLM::glmNet
#' @inheritParams glmnet_measure
#' @inheritParams measure_glm_raw
#' @param print_check Logical. When \code{TRUE}, prints intermediate results.
#' @return A list or a data frame. When \code{output_param_est = FALSE},
#' returns a data frame with a single row containing measures of model fitness.
#' Otherwise, returns a list with 2 elements. The first element,
#' \code{model_fitness}, contains a data frame with a single row containing
#' measures of model fitness, and the second element, \code{param_est},
#' contains a data frame of parameter estimates.
#' @examples
#' xtr <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' xte <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' b <- rnorm(5)
#'
#' ## continuous
#' ytr <- xtr %*% b + rnorm(100)
#' yte <- xte %*% b + rnorm(100)
#'
#' validate_ssnet(
#'   model = "glmnet", family = "gaussian",
#'   x.train = xtr, x.test = xte,
#'   y.train = ytr, y.test = yte
#' )
#'
#' validate_ssnet(
#'   model = "ss", family = "gaussian",
#'   x.train = xtr, x.test = xte,
#'   y.train = ytr, y.test = yte
#' )
#'
#'  ## binary
#' ybtr <- ifelse(ytr > 0, 1, 0)
#' ybte <- ifelse(yte > 0, 1, 0)
#'
#' validate_ssnet(
#'   model = "glmnet", family = "binomial",
#'   x.train = xtr, x.test = xte,
#'   y.train = ybtr, y.test = ybte,
#'   classify = TRUE, s0 = 0.1
#' )
#'
#' validate_ssnet(
#'   model = "ss", family = "binomial",
#'   x.train = xtr, x.test = xte,
#'   y.train = ybtr, y.test = ybte,
#'   classify = TRUE, s0 = 0.05, s1 = 1
#' )
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
#' validate_ssnet(
#'   model = "glmnet", family = "multinomial",
#'   x.train = xtr, x.test = xte,
#'   y.train = ymtr, y.test = ymte,
#'   classify = TRUE, s0 = 0.1,
#'   output_param_est = TRUE
#' )
#'
#' validate_ssnet(
#'   model = "ss", family = "multinomial",
#'   x.train = xtr, x.test = xte,
#'   y.train = ymtr, y.test = ymte,
#'   classify = TRUE, s0 = 0.1, s1 = 1,
#'   output_param_est = TRUE
#' )
#'
#'
#' @export
validate_ssnet <- function(
  model = "ss",
  alpha = 1,
  classify = FALSE,
  classify.rule = 0.5,
  type.multinomial = "grouped",
  s0 = 0.01,
  s1 = 1,
  x.train,
  y.train,
  x.test,
  y.test,
  family = "gaussian",
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
  adjmat = NULL,
  p.bound = c(0.01, 0.99),
  tau.prior = "none",
  tau.manual = NULL,
  stan_manual = NULL,
  nlambda = 100,
  lambda.criteria = "lambda.min",
  output_param_est = FALSE,
  output_probs = FALSE,
  print_check = FALSE
){
  if (!(model %in% c("glmnet", "ss", "ss_iar"))) {
    stop("Models must be one of glmnet, ss, ss_iar.")
  }
  if (!all(family %in% c("gaussian", "binomial", "poisson", "multinomial"))) {
    stop("Models must be a combination of gaussian, binomial, poisson, or multinomial.")
  }
  if (family == "binomial") {
      model_fit <- c("deviance", "mse", "mae", "auc", "misclassification")
      dispersion <- 1
    }
  if (family %in% c("gaussian", "poisson")){
      model_fit <- c("deviance", "mse", "mae")
      if (family == "poisson") {
        dispersion <- 1
      }
  }
  if (ncol(x.train) != ncol(x.test)) {
    stop("x.train and x.test must have the same number of columns.")
  }

  var_names <- c()
  for (xn in seq_len((ncol(x.train) + 1))) {
    var_names[xn] <- paste0("x", xn - 1)
  }

  if (model == "glmnet") {
    fit.mod <- glmnet::glmnet(
      x = x.train, y = y.train,
      family = family, alpha = alpha,
      type.multinomial = type.multinomial,
      lambda = s0
    )
    if (family == "gaussian") {
      dispersion <- NULL
      warning("dispersion for gaussian model not implemented (yet)")
    }
    if (family == "multinomial") {
      pred.y <- predict(
        object = fit.mod,
        newx = x.test,
        s = s0,
        type = "class"
      )
      b.list <- fit.mod$beta
      b.list2 <- list()
      for (i in seq_len(length(b.list))) {
        b.list2[[i]] <- as.numeric(b.list[[i]])
        names(b.list2[[i]]) <- var_names[-1]
      }
      names(b.list2) <- names(b.list)
      b <- b.list2
      # print(b)
      pred.pc <- prob_multinomial(
        x = x.test, b = b, a0 = fit.mod$a0
      )
      if (output_param_est == TRUE) {
        param_est <- NULL
        coef.fit.mod <- coef(fit.mod, s = s0)
        for (i in seq_len(length(unique(y.train)))){
          param_est_i <- data.frame(matrix(
            coef.fit.mod[[i]], nrow = 1
          ))
          colnames(param_est_i) <- var_names
          param_est <- dplyr::bind_rows(
            param_est,
            param_est_i
          )
        }
        param_est <- cbind(
          outcome = names(coef.fit.mod),
          param_est
        )
      }
    } else {
      pred.y.raw <- predict(
        object = fit.mod,
        newx = x.test,
        s = s0,
        type = "link"
      )
      if (output_param_est == TRUE) {
        param_est <- data.frame(matrix(coef(fit.mod, s = s0), nrow = 1))
        colnames(param_est) <- var_names
      }
    }
  }

  if (model %in% c("ss", "ss_iar")) {
    fit.mod <- ssnet(
      x = x.train,
      y = y.train,
      family = family,
      alpha = alpha,
      ss = c(s0, s1),
      type.multinomial = type.multinomial,
      iar.prior = iar.prior,
      adjmat = adjmat,
      iar.data = iar.data,
      tau.prior = tau.prior,
      tau.manual = tau.manual,
      stan_manual = stan_manual,
      opt.algorithm = opt.algorithm,
      p.bound = p.bound
    )
    if (family == "multinomial") {
      pred.pc <- prob_multinomial(
        x = x.test, b = fit.mod$coefficients, a0 = fit.mod$a0
      )
      pred.y <- pred.pc$predicted.class

      if (output_param_est == TRUE) {
        param_est <- NULL
        coef.fit.mod <- fit.mod$coefficients
        for (i in seq_len(length(unique(y.train)))){
          param_est_i <- c(
            x0 = fit.mod$a0[i],
            fit.mod$coefficients[[i]]
            )
          names(param_est_i) <- var_names
          param_est <- dplyr::bind_rows(
            param_est,
            param_est_i
          )
        }
        param_est <- cbind(
          outcome = names(coef.fit.mod),
          param_est
        )
      }
    } else {
      pred.y.raw <- cbind(1, x.test) %*% fit.mod$coefficients
      if (output_param_est == TRUE) {
        param_est <- data.frame(matrix(
          fit.mod$coefficients, nrow = 1
          ))
        colnames(param_est) <- var_names
      }
    }
  }
  if (family == "multinomial") {
    model_fitness <- measures_multiclass(
      y = y.test,
      pr_yi = subset(pred.pc, select = -predicted.class),
      # pr_yi = pred.pc %>% dplyr::select(-predicted.class),
      print_check = print_check
    )
  } else {
    if (family == "gaussian") {
      pred.y <- pred.y.raw
      dispersion <- fit.mod$dispersion
    }
    if (family == "binomial") {
      exp.pred.y.raw <- exp(pred.y.raw)
      pred.y <- exp.pred.y.raw / (1 + exp.pred.y.raw)
      dispersion <- 1
    }
    if (family == "poisson") {
      pred.y <- exp(pred.y.raw)
      dispersion <- 1
    }
    model_fitness <- measures_model_fitness(
      y = y.test,
      y.fitted = pred.y,
      family = family,
      dispersion = dispersion,
      inverse.link.y = TRUE,
      classify = classify,
      classify.rule = classify.rule
    )
  }
  if (output_param_est == FALSE & output_probs == FALSE) {
    return(model_fitness)
  } else {

  }
  if (output_param_est == TRUE) {
    out <- list(
      model_fitness = model_fitness,
      param_est = param_est
    )
    if (output_probs == TRUE & family == "multinomial") {
      out$cat_probs <- pred.pc
    }
    return(out)
  }
}
