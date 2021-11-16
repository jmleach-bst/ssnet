library(sim2Dpredictr)
library(ssnet)

set.seed(9038908)

block_correlation_builder <- function(
  size_block = 20,
  num_block = 50,
  corr_block = 0.90
) {
  cm <- list()
  for (i in 1:num_block) {
    cmi <- matrix(corr_block, nrow = size_block, ncol = size_block)
    diag(cmi) <- 1
    cm[[i]] <- cmi
  }
  bm <- Matrix::bdiag(cm)
  return(bm)
}

R <- block_correlation_builder(
  size_block = 5,
  num_block = 5,
  corr_block = 0.50
)

# if you want non-unit covariance ... probably not, but just in case
sigma <- 1
if (sigma == 1) {
  S <- R
} else {
  D <- diag(sigma, nrow = p, ncol = p)
  S <- D %*% R %*% D
}

# spam is better
S <- spam::as.spam(as.matrix(S))
Rc <- spam::chol.spam(S)
L <- spam::t(Rc)

BB <- rep(0, 25)

# non-zero locations
nzl <- c(1, 6, 11, 16, 21)

# "reference" group 1 (all 0's)
B1 <- c(0, BB)
# group 2
B2 <- BB
B02 <- -0.15
B2[nzl] <- 0.125
B2 <- c(B02, B2)
# group 3
B3 <- BB
B03 <- -0.5
B3[nzl] <- 0.25
B3 <- c(B03, B3)

B.all <- list(B1 = B1,
              B2 = B2,
              B3 = B3)

df.ex <- sim2Dpredictr::sim_Y_MVN_X(
  N = 100, B = B.all, L = L, S = S,
  dist = "multinomial", V = 3, incl.subjectID = FALSE
)

fit.ex <- ssnet(
  family = "multinomial", type.multinomial = "grouped",
  alpha = 1, ss = c(0.2, 5),
  x = as.matrix(df.ex |> dplyr::select(-Y)),
  y = df.ex$Y, print.iter = FALSE)
