#' Fit Several Models and Compare
#'
#' Fit \code{\link[glmnet]{glmnet}} and/or \code{\link[ssnet]{ssnet}} models and output measures of model fit for each. Allows multiple
#' scale values for spike-and-slab models.
#'
#' @importFrom glmnet glmnet
#' @param models A vector that determines which models to fit. Options include \code{c("glmnet", "ss", "ss_iar")}.
#' The default is to fit all three models.
#' @param s0 A vector of user-selected possible values for the spike scale parameter. The default
#' is \code{s0 = seq(0.01, 0.1, 0.01)}.
#' @param s1 A scalar defining the slab scale parameter.
#' @param output_param_est Logical. When \code{TRUE} adds an element to the output list that includes parameter estimates
#' for each model fit. Defaults to \code{FALSE}.
#' @param model_fit A vector containing measures of model fit to output. Options include  \code{c("deviance", "mse", "mae")}
#' for all models, and when \code{family = "binomial"}, also \code{c("auc", "misclassification")}. When \code{model_fit = "all"},
#' then all appropriate measures of model fit are output.
#' @param variable_selection Logical. When \code{TRUE}, outputs the false dicovery proportion (FDP), family-wise error (FWE),
#' and power for the model. Requires that parameter vector \code{B} be specified. Default is \code{FALSE}, and is only
#' appropriate for simulated data, when the true and false positives can be known.
#' @param type_error Determines whether models are selected based on training error (\code{"training"})
#' or k-fold cross validated estimates of prediction error (\code{"kcv"}). Defaults to \code{"kcv"}, which
#' is recommended because training error tends to underestimate the generalization error. See, e.g., Ch. 7 in
#' \insertCite{Hastie:2009}{ssnet}.
#' @param B When \code{variable_selection} is \code{TRUE}, a vector of "true" parameter values must be input in order
#' to calculate the false discovery proportion (FDR), family-wise error (FWER), and power for the data. This vector should
#' NOT contain a value for the intercept.
#' @inheritParams ssnet
#' @inheritParams glmnet::glmnet
#' @inheritParams glmnet::cv.glmnet
#' @inheritParams BhGLM::glmNet
#' @inheritParams glmnet_measure
#' @inheritParams measure_glm_raw
#' @param criteria Specifies the criteria for model selection. Options are \code{"deviance"}, \code{"mse"},
#' \code{"mae"} for deviance, mean-square error, and mean absolute error, respectively. When
#' \code{family = "binomial"}, additional options are \code{"auc"} and \code{"misclassification"}, for
#' Area under the ROC curve and the percentage of cases where the difference between the observed and
#' predicted values is greater than 1/2.
#' @return When \code{output_param_est = FALSE} returns a data frame of model fitness summaries. Otherwise, returns
#' a list whose first element is a dataframe whose rows contain parameter estimates for each model fit, and whose
#' second element is a dataframe of model fitness summaries.
#' @note Models fit with `glmnet` never select the penalty/tuning parameter using the training error; however, when
#' \code{type_error = "training"}, the measure used to compare `glmnet` with the other models is based on prediction
#' error estimates from training error. That is, model selection within `glmnet` is still based on k-fold cross validation,
#' even if comparisons with other models is not.
#' @examples
#' library(sim2Dpredictr)
#' library(ssnet)
#' ## generate data (no intercept)
#' set.seed(4799623)
#' cn <- c()
#' for (i in 1:100) cn[i] <- paste0("x", i)
#' tb <- rbinom(100, 1, 0.05)
#' tx <- matrix(rnorm(10000), nrow = 100, ncol = 100,
#'              dimnames = list(1:100, cn))
#' ty <- tx %*% tb + rnorm(100)
#'
#' ## build adjacency matrix
#' adjmat <- proximity_builder(im.res = c(10, 10), type = "sparse")
#' ## stan model information
#' model_info <- mungeCARdata4stan(adjmat$nb.index,
#'                                 table(adjmat$location.index))
#' ## pre-specify stan model
#' sm <- stan_model(file =
#' "C:/Users/Justin/Documents/BST/Dissertation_in_Latex/stan models/iar_incl_prob_notau.stan"
#' )
#'
#' ## fit multiple models and compare
#' compare_ssnet(x = tx, y = ty, alpha = c(0, 0.5, 1),
#'               s0 = c(0.01, 0.05, 0.1),
#'               family = "gaussian", type_error = "kcv",
#'               model_fit = "all", variable_selection = TRUE,
#'               B = tb, iar.data = model_info, stan_manual = sm)
#' @export
compare_ssnet <- function(models = c("glmnet", "ss", "ss_iar"),
                          alpha = c(0, 0.5, 1),
                          model_fit = "all", variable_selection = FALSE,
                          classify = FALSE, classify.rule = 0.5,
                          type_error = "kcv", nfolds = 10, ncv = 1,
                          s0 = seq(0.01, 0.2, 0.02), s1 = 1, B = NULL,
                          x, y, family = "gaussian",
                          offset = NULL, epsilon = 1e-04, maxit = 50, init = NULL, group = NULL,
                          Warning = FALSE, verbose = FALSE, opt.algorithm = "LBFGS",
                          iar.data = NULL, iar.prior = FALSE,
                          p.bound = c(0.01, 0.99), tau.prior = "none",
                          stan_manual = NULL, stan_local = FALSE,
                          plot.pj = FALSE, im.res = NULL,
                          nlambda = 100,
                          type.measure = c("default", "mse", "deviance",
                                           "class", "auc", "mae", "C"),
                          lambda.criteria = "lambda.min",
                          output_param_est = FALSE) {

  options(stringsAsFactors = FALSE)
  if (variable_selection == TRUE) {
    if (is.null(B) == TRUE | length(B) != ncol(x)) {
      stop("Variable selection measures requires a true parameter vector of appropriate length.")
    }
  }
  if (!all(models %in% c("glmnet", "ss", "ss_iar"))) {
    stop("Models must be a combination of all, glmnet, ss, ss_iar. \n")
  }
  if (!all(family %in% c("gaussian", "binomial", "poisson", "cox"))) {
    stop("Models must be a combination of gaussian, binomial, poisson, cox. \n")
  }
  if (model_fit == "all") {
    if (family == "binomial") {
      model_fit = c("deviance", "mse", "mae", "auc", "misclassification")
    } else {
      model_fit = c("deviance", "mse", "mae")
    }
  }

  out <- NULL
  glmnet.fit.a <- NULL

  # variable names
  xi <- c("x0")
  for (i in 1:ncol(x)) {
    xi <- c(xi, paste0("x", i))
  }

  # fit glmnet
  if ("glmnet" %in% models) {

        glmnet.params <- NULL

    for (a in 1:length(alpha)) {
      if (type_error == "training") {
        glmnet.fit.a <- glmnet::cv.glmnet(x = x, y = y, family = family, alpha = alpha[a],
                                          type.measure = type.measure, nlambda = nlambda,
                                          nfolds = nfolds)
        glmnet.mp.a <- glmnet_measure(glmnet.fit.a, x = x, y = y, family = family,
                                      type.measure = type.measure,
                                      lambda.criteria = lambda.criteria, alpha = alpha[a])
        if (output_param_est == TRUE) {
          glmnet.params.a <- as.vector(glmnet::coef.glmnet(glmnet.fit.a, s = lambda.criteria))
          glmnet.params.a <- data.frame(model = "glmnet",
                                      alpha = alpha[a],
                                      s0 = glmnet.fit.a[[lambda.criteria]],
                                      t(glmnet.params.a))
          names(glmnet.params.a) <- c("model", "alpha", "s0", xi)
          if (a == 1) {
            glmnet.params <- glmnet.params.a
          } else {
            glmnet.params <- rbind(glmnet.params,
                                   glmnet.params.a)
          }
          # print(glmnet.params)
        }
      } else {
        glmnet.fit.a <- BhGLM::glmNet(x = x, y = y, family = family,
                                    alpha = alpha[a],
                                    nfolds = nfolds, ncv = ncv,
                                    verbose = verbose)
        if (ncv == 1) {
          glmnet.mp.a <- as.data.frame(t(cv.bh3(glmnet.fit.a, nfolds = nfolds, ncv = ncv,
                                                classify = classify, classify.rule = classify.rule,
                                                verbose = verbose)$measures))
        } else {
          glmnet.mp.a <- as.data.frame(t(cv.bh3(glmnet.fit.a, nfolds = nfolds, ncv = ncv,
                                                classify = classify, classify.rule = classify.rule,
                                                verbose = verbose)$measures[1, ]))
        }
        glmnet.mp.a <- cbind(model = "glmnet", alpha = alpha[a], s0 = glmnet.fit.a$lambda, glmnet.mp.a)
        if (output_param_est == TRUE) {
          glmnet.params.a <- data.frame(model = "glmnet",
                                        alpha = alpha[a],
                                        s0 = glmnet.fit.a$lambda,
                                        t(glmnet.fit.a$coefficients))
          names(glmnet.params.a) <- c("model", "alpha", "s0", xi)
          glmnet.params <- rbind(glmnet.params, glmnet.params.a)
        }
      }

      if (variable_selection == TRUE) {
        glmnet.rej.a <- rep(0, ncol(x))
        glmnet.est.a <- stats::coef(glmnet.fit.a, s = lambda.criteria)[-1]
        glmnet.rej.a[glmnet.est.a != 0] = 1
        glmnet.vs.a <- sim2Dpredictr::sample_FP_Power(rejections = glmnet.rej.a, B = B, B.incl.B0 = FALSE)
        glmnet.mp.a <- cbind(glmnet.mp.a, glmnet.vs.a)
      }
      if (a == 1) {
        glmnet.mp <- glmnet.mp.a
      } else {
        glmnet.mp <- rbind(glmnet.mp,
                           glmnet.mp.a)
      }
      # print(glmnet.mp)
    }
  }

  if (is.null(glmnet.fit.a) == FALSE) {
    out <- glmnet.mp
    if (output_param_est == TRUE) {
      param.est <- glmnet.params
    }
  } else {
    out <- NULL
    if (output_param_est == TRUE) {
      param.est <- NULL
    }
  }

  # fit ssnet
  if ( "ss" %in% models) {

    ss.mp <- NULL
    ss.params <- NULL
    ss.fit.i.a <- NULL
    ssiar.fit.i.a <- NULL

    for (i in 1:length(s0)) {
      ss <- c(s0[i], s1)
      for (a in 1:length(alpha)) {
        ss.fit.i.a <- ssnet::ssnet(x = x, y = y, family = family, ss = ss,
                                        offset = offset, init = init,
                                        group = group, alpha = alpha[a],
                                        iar.prior = FALSE,
                                        verbose = verbose)
        if(output_param_est == TRUE) {
          ss.params.a <- data.frame(model = "ss",
                                    alpha = alpha[a],
                                    s0 = s0[i],
                                    t(ss.fit.i.a$coefficients))
          ss.params <- rbind(ss.params, ss.params.a)
        }

        # model fit measures
        if (type_error == "training") {
          ss.mp.i.a <- as.data.frame(t(c(alpha = alpha[a], s0 = s0[i], measure_bh_raw(ss.fit.i.a)[model_fit])))
        } else {
          if (ncv == 1) {
            ss.mp.i.a <- as.data.frame(t(c(alpha = alpha[a],
                                           s0 = s0[i],
                                           cv.bh3(ss.fit.i.a,
                                                  nfolds = nfolds,
                                                  ncv = ncv,
                                                  classify = classify,
                                                  classify.rule = classify.rule,
                                                  verbose = verbose)$measures)))
          } else {
            ss.mp.i.a <- as.data.frame(t(c(alpha = alpha[a],
                                           s0 = s0[i],
                                           cv.bh3(ss.fit.i.a,
                                                  nfolds = nfolds,
                                                  ncv = ncv,
                                                  classify = classify,
                                                  classify.rule = classify.rule,
                                                  verbose = verbose)$measures[1, ])))
          }
        }
        # varaible selection performance
        if (variable_selection == TRUE) {
          ss.rej.i.a <- rep(0, ncol(x))
          ss.est.i.a = ss.fit.i.a$coefficients[-1]
          ss.rej.i.a[ss.est.i.a != 0] = 1
          ss.vs.i.a <- sim2Dpredictr::sample_FP_Power(rejections = ss.rej.i.a, B = B, B.incl.B0 = FALSE)
          ss.mp.i.a <- cbind(ss.mp.i.a, ss.vs.i.a)
        }
        ss.mp <- rbind(ss.mp, ss.mp.i.a)
      }
    }

    ss.mp <- data.frame(model = rep("ss", nrow(ss.mp)), ss.mp)
    #ss.mp$model <- "ss"
    row.names(ss.mp) <- NULL

    # combine with previous models
    # print(ss.mp)
    if (is.null(out) == FALSE) {
      out <- rbind(out, ss.mp)
      if (output_param_est == TRUE) {
        names(ss.params) <- c("model", "alpha", "s0", xi)
        param.est <- rbind(glmnet.params, ss.params)
      }
    } else {
      out <- ss.mp
      if (output_param_est == TRUE) {
        names(ss.params) <- c("model", "alpha", "s0", xi)
        param.est <- ss.params
      }
    }
  }

  # fit ssnet_iar
  if ( "ss_iar" %in% models) {

    ssiar.mp <- NULL
    ssiar.params <- NULL

    for (i in 1:length(s0)) {
      ss <- c(s0[i], s1)
      for (a in 1:length(alpha)) {
        ssiar.fit.i.a <- ssnet::ssnet(x = x, y = y, family = family, ss = ss,
                                        offset = offset, init = init,
                                        group = group, alpha = alpha[a],
                                        iar.prior = TRUE,
                                        verbose = verbose,
                                        iar.data = iar.data, p.bound = p.bound,
                                        tau.prior = tau.prior,
                                        stan_manual = stan_manual,
                                        stan_local = stan_local,
                                        im.res = im.res)
        if(output_param_est == TRUE) {
          ssiar.params.a <- data.frame(model = "ss_iar",
                                       alpha = alpha[a],
                                       s0 = s0[i],
                                       t(ssiar.fit.i.a$coefficients))
          ssiar.params <- rbind(ssiar.params, ssiar.params.a)
        }

        # model fit measures
        if (type_error == "training") {
          ssiar.mp.i.a <- as.data.frame(t(c(alpha = alpha[a], s0 = s0[i], measure_bh_raw(ssiar.fit.i.a)[model_fit])))
        } else {
          if (ncv == 1) {
            ssiar.mp.i.a <- as.data.frame(t(c(alpha = alpha[a],
                                           s0 = s0[i],
                                           cv.bh3(ssiar.fit.i.a,
                                                  nfolds = nfolds,
                                                  ncv = ncv,
                                                  classify = classify,
                                                  classify.rule = classify.rule,
                                                  verbose = verbose)$measures)))
          } else {
            ssiar.mp.i.a <- as.data.frame(t(c(alpha = alpha[a],
                                           s0 = s0[i],
                                           cv.bh3(ssiar.fit.i.a,
                                                  nfolds = nfolds,
                                                  ncv = ncv,
                                                  classify = classify,
                                                  classify.rule = classify.rule,
                                                  verbose = verbose)$measures[1, ])))
          }
        }
        # varaible selection performance
        if (variable_selection == TRUE) {
          ssiar.rej.i.a <- rep(0, ncol(x))
          ssiar.est.i.a = ssiar.fit.i.a$coefficients[-1]
          ssiar.rej.i.a[ssiar.est.i.a != 0] = 1
          ssiar.vs.i.a <- sim2Dpredictr::sample_FP_Power(rejections = ssiar.rej.i.a, B = B, B.incl.B0 = FALSE)
          ssiar.mp.i.a <- cbind(ssiar.mp.i.a, ssiar.vs.i.a)
        }
        ssiar.mp <- rbind(ssiar.mp, ssiar.mp.i.a)
      }
    }

    ssiar.mp <- data.frame(model = rep("ss_iar", nrow(ssiar.mp)), ssiar.mp)
    #ssiar.mp$model <- "ss"
    row.names(ssiar.mp) <- NULL

    # combine with previous models
    # print(ssiar.mp)
    if (is.null(out) == FALSE) {
      out <- rbind(out, ssiar.mp)
      if (output_param_est == TRUE) {
        names(ssiar.params) <- c("model", "alpha", "s0", xi)
        param.est <- rbind(param.est, ssiar.params)
      }
    } else {
      out <- ssiar.mp
      if (output_param_est == TRUE) {
        names(ssiar.params) <- c("model", "alpha", "s0", xi)
        param.est <- ssiar.params
      }
    }
  }

  if (output_param_est == TRUE) {
    return(list(estimates = param.est, inference = out))
  } else {
    return(out)
  }
}
