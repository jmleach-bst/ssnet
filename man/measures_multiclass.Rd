% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_multiclass.R
\name{measures_multiclass}
\alias{measures_multiclass}
\title{Obtain Measures of Model Performance for Mutliclass Classification}
\usage{
measures_multiclass(y, pr_yi = NULL, y_hat = NULL, print_check = FALSE)
}
\arguments{
\item{y}{Scalar response variable denoting class membership. It is preferred
that \code{y} is a factor with two or more levels, but will modify
internally for numeric or character variables.}

\item{pr_yi}{A data frame or matrix where each row is a subject observation
where columns indicate the probability the subjects belongs to the
respective class. The column names should match the factor levels for
\code{y}. This information can be obtained using \code{prob_multinomial()}.}

\item{y_hat}{A vector of predicted classes that should have the same factor
levels as \code{y}. When \code{y_hat} is supplied it supersedes
\code{pr_yi} for classification. However, \code{pr_yi} will still be used
to calculate the deviance.}

\item{print_check}{Logical. When \code{TRUE}, prints intermediate results.}
}
\value{
A data frame with a single row containing columns for average
accuracy (avg_acc), average per-class classification error (pce),
micro-averaged positive predictive value (ppv_micro), sensitivity (sn_micro),
F1 score (f1_micro), and macro-averaged positive predictive value (ppv_macro),
sensitivity (sn_macro), F1 score (f1_macro).
}
\description{
Compare observed and predicted classes in cases where there are more than
2 classes.
}
\details{
Multiclass measures are calculated according to Sokolova and
Lapalme (2009). We give brief details here, but refer users to Table 3 of
Sokolova & Lapalme (2009) for precise definitions.
}
\note{
While this function will work when there are only 2 classes, it is
recommended to use \code{eval_classify()} when there are only 2 classes.
}
\examples{

n <- 500
y <- sample(c("ack", "eek", "ahh"), size = n, replace = TRUE)
x <- matrix(rnorm(n*5), nrow = n, ncol = 5)
colnames(x) <- paste0("x", seq_len(ncol(x)))
mfit <- glmnet::glmnet(x = x, y = y, family = "multinomial",
                       type.multinomial = "grouped",
                       lambda = 0.01)
xnew <- matrix(rnorm(n*5), nrow = n, ncol = 5)
b <- list()
for (i in seq_len(length(mfit$beta))) {
  b[[i]] <- as.numeric(mfit$beta[[i]])
}
names(b) <- names(mfit$beta)
for (i in seq_len(length(b))) {
    names(b[[i]]) <- colnames(x)
}
pm <- prob_multinomial(x = xnew, b = b, a0 = mfit$a0)
measures_multiclass(y = y, pr_yi = pm[,-4])
measures_multiclass(y = y, y_hat = pm$predicted.class)

}
\references{
\insertRef{Sokolova+Lapalme:2009}{ssnet}
}
