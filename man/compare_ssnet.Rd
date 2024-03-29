% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_ssnet.R
\name{compare_ssnet}
\alias{compare_ssnet}
\title{Fit Several Models and Compare}
\usage{
compare_ssnet(
  models = c("glmnet", "ss", "ss_iar"),
  alpha = c(0, 0.5, 1),
  model_fit = "all",
  variable_selection = FALSE,
  classify = FALSE,
  classify.rule = 0.5,
  type_error = "kcv",
  nfolds = 10,
  ncv = 1,
  foldid = NULL,
  fold.seed = NULL,
  s0 = seq(0.01, 0.1, 0.01),
  s1 = 1,
  B = NULL,
  x,
  y,
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
  p.bound = c(0.01, 0.99),
  tau.prior = "none",
  stan_manual = NULL,
  plot.pj = FALSE,
  im.res = NULL,
  nlambda = 100,
  type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
  lambda.criteria = "lambda.min",
  output_param_est = FALSE,
  output_cv = FALSE
)
}
\arguments{
\item{models}{A vector that determines which models to fit. Options include
\code{c("glmnet", "ss", "ss_iar")}. The default is to fit all three models.}

\item{alpha}{A scalar value between 0 and 1 determining the compromise
between the Ridge and Lasso models. When \code{alpha = 1} reduces to the
Lasso, and when \code{alpha = 0} reduces to Ridge.}

\item{model_fit}{A vector containing measures of model fit to output.
Options include  \code{c("deviance", "mse", "mae")} for all models, and
when \code{family = "binomial"}, also \code{c("auc", "misclassification")}.
When \code{model_fit = "all"}, then all appropriate measures of model fit
are output.}

\item{variable_selection}{Logical. When \code{TRUE}, outputs the false
discovery proportion (FDP), family-wise error (FWE), and power for the
model. Requires that parameter vector \code{B} be specified. Default is
\code{FALSE}, and is only appropriate for simulated data, when the true and
false positives can be known.}

\item{classify}{Logical. When \code{TRUE} and \code{family = "binomial"}
applies a classification rule given by the argument \code{classify.rule},
and outputs accuracy, sensitivity, specificity, positive predictive value
(ppv), and negative predictive value (npv).}

\item{classify.rule}{A value between 0 and 1. For a given predicted value
from a logistic regression, if the value is above \code{classify.rule},
then the predicted class is 1; otherwise the predicted class is 0. The
default is 0.5.}

\item{type_error}{Determines whether models are selected based on training
error (\code{"training"}) or k-fold cross validated estimates of prediction
error (\code{"kcv"}). Defaults to \code{"kcv"}, which is recommended because
training error tends to underestimate the generalization error. See, e.g.,
Ch. 7 in \insertCite{Hastie:2009}{ssnet}.}

\item{nfolds}{number of folds - default is 10. Although \code{nfolds} can be
as large as the sample size (leave-one-out CV), it is not recommended for
large datasets. Smallest value allowable is \code{nfolds=3}}

\item{ncv}{
  repeated number of cross-validation.  
  }

\item{foldid}{an optional vector of values between 1 and \code{nfolds}
identifying what fold each observation is in. If supplied, \code{nfolds} can
be missing.}

\item{fold.seed}{A scalar seed value for cross validation; ensures the
folds are the same upon re-running the function. Alternatively, use
\code{foldid} to manually specify folds.}

\item{s0, s1}{A vector of user-selected possible values for the spike scale
and  slab scale parameter, respectively. The default is
\code{s0 = seq(0.01, 0.1, 0.01)} and \code{s1 = 1}. However, the user
should select values informed by the practical context of the analysis.}

\item{B}{When \code{variable_selection} is \code{TRUE}, a vector of "true"
parameter values must be input in order to calculate the false discovery
proportion (FDR), family-wise error (FWER), and power for the data. This
vector should NOT contain a value for the intercept.}

\item{x}{Design, or input, matrix, of dimension nobs x nvars; each row is
an observation vector. It is recommended that \code{x} have user-defined
column names for ease of identifying variables. If missing, then
\code{colnames} are internally assigned \code{x1}, \code{x2}, ... and so
forth.}

\item{y}{Scalar response variable. Quantitative for
\code{family = "gaussian"}, or \code{family = "poisson"}
(non-negative counts). For \code{family = "gaussian"}, \code{y} is always
standardized. For \code{family = "binomial"}, \code{y} should be either a
factor with two levels, or a two-column matrix of counts or proportions
(the second column is treated as the target class; for a factor, the last
level in alphabetical order is the target class). For \code{family="cox"},
\code{y} should be a two-column matrix with columns named \code{'time'}
and \code{'status'}. The latter is a binary variable, with \code{'1'}
indicating death, and \code{'0'} indicating right censored. The function
\code{Surv()} in package survival produces such a matrix. When
\code{family = "multinomial"}, \code{y} follows the documentation for
\code{glmnet}, but it is preferred that \code{y} is a factor with two
or more levels.}

\item{family}{Response type (see above).}

\item{offset}{A vector of length \code{nobs} that is included in the linear
predictor.}

\item{epsilon}{A positive convergence tolerance; the iterations converge
when \eqn{|dev - dev_old|/(|dev| + 0.1) < e}.}

\item{maxit}{An integer giving the maximal number of EM iterations.}

\item{init}{A vector of initial values for all coefficients (not for
intercept). If not given, it will be internally produced. If
\code{family = "multinomial"} and the same initializations are desired for
each response/outcome category then \code{init} can be a vector. If
different initializations are desired, then \code{init} should be a list,
each element of which contains a vector of initializations. The list should
be named according the response/outcome category as they appear in \code{y}.}

\item{group}{A numeric vector, or an integer, or a list indicating the
groups of predictors. If \code{group = NULL}, all the predictors form a
single group. If \code{group = K}, the predictors are evenly divided into
groups each with K predictors. If group is a numberic vector, it defines
groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2:
\code{(group[2]+1):group[3]}, Group 3: \code{(group[3]+1):group[4]}, ...
If group is a list of variable names, \code{group[[k]]} includes variables
in the k-th group. The mixture double-exponential prior is only used for
grouped predictors. For ungrouped predictors, the prior is
double-exponential with scale \code{ss[2]} and mean 0. Note that grouped
predictors when \code{family = "multinomial"} is still experimental, so
use with caution.}

\item{Warning}{Logical. If \code{TRUE}, shows the error messages of not
convergence and identifiability.}

\item{verbose}{Logical. If \code{TRUE}, prints out the number of iterations
and computational time.}

\item{opt.algorithm}{One of \code{c("LBFGS", "BFGS", "Newton")}. This
argument determines which argument is used to optimize the term in the EM
algorithm that estimates the probabilities of inclusion for each parameter.
Optimization is performed by \code{optimizing}.}

\item{iar.data}{A list of output from \code{\link{mungeCARdata4stan}} that
contains the necessary inputs for the IAR prior. When unspecified, this is
built internally assuming that neighbors are those variables directly above,
below, left, and  right of a given variable location. \code{im.res} must be
specified when allowing this argument to be built internally. It is not
recommended to use this argument directly, even when specifying a more
complicated neighborhood stucture; this can be specified with the
\code{adjmat} argument, and then internally converted to the correct format.}

\item{iar.prior}{Logical. When \code{TRUE}, imposes intrinsic autoregressive
prior on logit of the probabilities of inclusion. When \code{FALSE}, treats
probabilities of inclusion as unstructured.}

\item{p.bound}{A vector defining the lower and upper boundaries for the
probabilities of inclusion in the model, respectively. Defaults to
\code{c(0.01, 0.99)}.}

\item{tau.prior}{One of \code{c("none", "manual", "cauchy")}. This argument
determines the precision parameter in the Conditional Autoregressive model
for the (logit of) prior inclusion probabilities. When \code{"none"}, the
precision is set to 1; when "manual", the precision is manually entered by
the user; when \code{"cauchy"}, the inverse precision is assumed to follow
a Cauchy distribution with mean 0 and scale 2.5. Note that at this stage of
development, only the \code{"none"} option has been extensively tested, so
the other options should be used with caution.}

\item{stan_manual}{A \code{stan_model} that is manually specified.
Especially when fitting multiple models in succession, specifying the
\code{stan} model outside this "loop" may avoid errors.}

\item{plot.pj}{When \code{TRUE}, prints a series of 2D graphs of the prior
probabilities of inclusion at each step of the algorithm. This should NOT
be used for 3D data.}

\item{im.res}{A 2-element vector where the first argument is the number of
"rows" and the second argument is the number of "columns" in each subject's
"image". Default is \code{NULL}.}

\item{nlambda}{The number of \code{lambda} values - default is 100.}

\item{type.measure}{loss to use for cross-validation. Currently five
options, not all available for all models. The default is
\code{type.measure="deviance"}, which uses squared-error for gaussian models
(a.k.a \code{type.measure="mse"} there), deviance for logistic and poisson
regression, and partial-likelihood for the Cox model.
\code{type.measure="class"} applies to binomial and multinomial logistic
regression only, and gives misclassification error.
\code{type.measure="auc"} is for two-class logistic regression only, and
gives area under the ROC curve. \code{type.measure="mse"} or
\code{type.measure="mae"} (mean absolute error) can be used by all models
except the \code{"cox"}; they measure the deviation from the fitted mean to
the response.
\code{type.measure="C"} is Harrel's concordance measure, only available for \code{cox} models.}

\item{lambda.criteria}{Determines the model selection criteria. When
\code{"lambda.min"} the final model is selected based on the penalty that
minimizes the measure given in \code{type.measure}. When \code{"lambda.1se"}
the final model is selected based on the smallest value of lambda that is
within one standard error of the minimal measure given in
\code{type.measure}.}

\item{output_param_est}{Logical. When \code{TRUE} adds an element to the
output list that includes parameter estimates for each model fit. Defaults
to \code{FALSE}.}

\item{output_cv}{Logical. When \code{TRUE}, generates a list with an
element for all output generated by \code{cv.bh3}. Default is \code{FALSE}.
The master list has an element for every combination of model, alpha, and
s0. Each element is itself a list, the first element being the model
fitness measures, and the second element the observed outcomes, fitted
outcomes, linear predictor, and fold identifier. When \code{ncv > 1}, the
latter three will consist of multiple columns, one for each run of k-fold
CV; e.g., the column containing the fitted outcomes for the second CV run
is \code{y.fitted_2}.}

\item{criteria}{Specifies the criteria for model selection. Options are
\code{"deviance"}, \code{"mse"}, \code{"mae"} for deviance, mean-square
error, and mean absolute error, respectively. When
\code{family = "binomial"}, additional options are \code{"auc"} and
\code{"misclassification"}, for area under the ROC curve and the percentage
of cases where the difference between the observed and predicted values is
greater than 1/2.}
}
\value{
When \code{output_param_est = FALSE} returns a data frame of model
fitness summaries. Otherwise, returns a list whose first element is a
dataframe whose rows contain parameter estimates for each model fit, and
whose second element is a dataframe of model fitness summaries.
}
\description{
Fit \code{\link[glmnet]{glmnet}} and/or \code{\link[ssnet]{ssnet}} models
and output measures of model fit for each. Allows multiple scale values for
spike-and-slab models.
}
\note{
Models fit with `glmnet` never select the penalty/tuning parameter
using the training error; however, when \code{type_error = "training"},
the measure used to compare `glmnet` with the other models is based on
prediction error estimates from training error. That is, model selection
within `glmnet` is still based on k-fold cross validation, even if
comparisons with other models is not.
}
\examples{
library(sim2Dpredictr)
## generate data (no intercept)
set.seed(4799623)

## sample size
n <- 30
## image dims
nr <- 4
nc <- 4

## generate data
cn <- paste0("x", seq_len(nr * nc))
tb <- rbinom(nr * nc, 1, 0.05)
tx <- matrix(rnorm(n * nr * nc), nrow = n, ncol = nr * nc,
             dimnames = list(seq_len(n), cn))
ty <- tx \%*\% tb + rnorm(n)

## build adjacency matrix
adjmat <- proximity_builder(im.res = c(nr, nc), type = "sparse")

## fit multiple models and compare
compare_ssnet(x = tx, y = ty, family = "gaussian",
              alpha = c(0.5, 1), s0 = c(0.01, 0.05),
              type_error = "kcv", nfolds = 3, im.res = c(4, 4),
              model_fit = "all", variable_selection = TRUE,
              B = tb)

}
