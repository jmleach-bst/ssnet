test_that("eval_classify throws correct errors", {
  expect_error(eval_classify(num.class = 1),
               "Require at least 2 classes.")
  expect_error(eval_classify(num.class = 3, beta.hat = c("a", "b")),
               "Parameter vector beta.hat must contain only numeric values.")
  expect_error(eval_classify(num.class = 3, beta.hat = 1:5, y = "1"),
               "y must contain only numeric values.")
  expect_error(eval_classify(num.class = 2, beta.hat = 1:3,
                             y = rnorm(10),
                             x = data.frame(x1 = rnorm(10),
                                            x2 = rnorm(10),
                                            x3 = sample(x = c("a", "b"),
                                                        size = 10,
                                                        replace = TRUE))),
               "x must contain only numeric values.")

})

test_that("vector and data frame inputs are equal", {
  expect_equal(
    eval_classify(
      beta.hat = data.frame(
        a = 1, b1 = 3, b2 = -3
        ),
      y = data.frame(
        s1 = 1, s2 = 0, s3 = 0,
        s4 = 1, s5 = 0, s6 = 1
        ),
      x = data.frame(x1 = 1:6, x2 = 1:6),
      family = "binomial",
      classify.rule = 0.5),
    eval_classify(
      beta.hat = c(a = 1, b1 = 3, b2 = -3),
      y = c(s1 = 1, s2 = 0, s3 = 0, s4 = 1, s5 = 0, s6 = 1),
      x = matrix(rep(1:6, 2), ncol = 2),
      family = "binomial",
      classify.rule = 0.5
      )
  )
})
