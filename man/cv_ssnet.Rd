% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_ssnet.R
\name{cv_ssnet}
\alias{cv_ssnet}
\title{Cross Validation for ssnet Models}
\usage{
cv_ssnet(
  model,
  alpha = c(0.5, 1),
  s0 = seq(0.01, 0.1, 0.01),
  s1 = c(1, 2.5),
  classify = FALSE,
  classify.rule = 0.5,
  nfolds = 10,
  ncv = 1,
  foldid = NULL,
  fold.seed = NULL,
  x,
  y,
  family,
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
  lambda.criteria = "lambda.min",
  output_param_est = FALSE,
  type.multinomial = "grouped"
)
}
\arguments{
\item{model}{Specify which model to fit. Options include
\code{c("glmnet", "ss", "ss_iar")}.}

\item{alpha}{A scalar value between 0 and 1 determining the compromise
between the Ridge and Lasso models. When \code{alpha = 1} reduces to the
Lasso, and when \code{alpha = 0} reduces to Ridge.}

\item{s0, s1}{A vector of user-selected possible values for the spike scale
and slab scale parameter, respectively. The default is
\code{s0 = seq(0.01, 0.1, 0.01)} and \code{s1 = 1}. However, the user
should select values informed by the practical context of the analysis.}

\item{classify}{Logical. When \code{TRUE} and \code{family = "binomial"}
applies a classification rule given by the argument \code{classify.rule},
and outputs accuracy, sensitivity, specificity, positive predictive value
(ppv), and negative predictive value (npv).}

\item{classify.rule}{A value between 0 and 1. For a given predicted value
from a logistic regression, if the value is above \code{classify.rule},
then the predicted class is 1; otherwise the predicted class is 0. The
default is 0.5.}

\item{nfolds}{Numeric value indicating the number of folds to create.}

\item{ncv}{Numeric value indicating the number of times to perform cross
validation.}

\item{foldid}{An (optional) vector of values between 1 and \code{nfold}
identifying the fold for each observation. When supplied \code{nfolds} may
be omitted. If \code{ncv > 1}, then supply a matrix or data frame  where
each column contains fold identifiers. If \code{foldid} is supplied, it
supersedes \code{ncv} and \code{nfolds}.}

\item{fold.seed}{A scalar seed value for cross validation; ensures the
folds are the same upon re-running the function. Alternatively, use
\code{foldid} to manually specify folds.}

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

\item{lambda.criteria}{Determines the model selection criteria. When
\code{"lambda.min"} the final model is selected based on the penalty that
minimizes the measure given in \code{type.measure}. When \code{"lambda.1se"}
the final model is selected based on the smallest value of lambda that is
within one standard error of the minimal measure given in
\code{type.measure}.}

\item{output_param_est}{Logical. When \code{TRUE} adds an element to the
output list that includes parameter estimates for the fitted model.
Defaults is \code{FALSE}.}

\item{type.multinomial}{If \code{"grouped"} then a grouped lasso penalty is
used on the multinomial coefficients for a variable. This ensures they are
all in our out together. The default is \code{"ungrouped"}}
}
\value{
Either a data frame of model fitness measures or a list whose
elements are data frames of model fitness measures and parameter estimates,
respectively, depending on the value of output_param_ets.
}
\description{
Perform k-fold cross validation for spike-and-slab elastic net models.
}
\examples{
xtr <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
xte <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
b <- rnorm(5)

## continuous outcome
ytr <- xtr \%*\% b + rnorm(100)
yte <- xte \%*\% b + rnorm(100)

## binary outcome
ybtr <- ifelse(ytr > 0, 1, 0)
ybte <- ifelse(yte > 0, 1, 0)

## multinomial outcome
ymtr <- dplyr::case_when(
  ytr > 1 ~ "a",
  ytr <= 1 & ytr > -1 ~ "b",
  ytr <= -1 ~ "c"
)
ymte <- dplyr::case_when(
  yte > 1 ~ "a",
  yte <= 1 & yte > -1 ~ "b",
  yte <= -1 ~ "c"
)

cv_ssnet(
  model = "ss", family = "gaussian",
  x = rbind(xtr, xte), y = c(ytr, yte),
  s0 = c(0.01, 0.05, 0.10), s1 = c(1, 2.5),
  nfolds = 3, ncv = 2
)

\dontrun{
cv_ssnet(
  model = "ss", family = "binomial",
  x = rbind(xtr, xte), y = c(ybtr, ybte),
  s0 = c(0.01, 0.05, 0.10), s1 = c(1, 2.5),
  nfolds = 3, ncv = 2, classify = TRUE,
  output_param_est = TRUE
)

cv_ssnet(
  model = "ss", family = "multinomial",
  x = rbind(xtr, xte), y = c(ymtr, ymte),
  s0 = c(0.01, 0.05, 0.10), s1 = c(1, 2.5),
  nfolds = 3, ncv = 2, classify = FALSE,
  output_param_est = TRUE
)
}

}
