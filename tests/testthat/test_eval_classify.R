test_that("vector and data frame inputs are equal", {
  expect_equal(eval_classify(beta.hat = data.frame(a = 1, b1 = 3, b2 = -3),
                             y = data.frame(s1 = 1, s2 = 0, s3 = 0, s4 = 1, s5 = 0, s6 = 1),
                             x = data.frame(x1 = 1:6, x2 = 1:6),
                             family = "binomial",
                             classify.rule = 0.5),
               eval_classify(beta.hat = c(a = 1, b1 = 3, b2 = -3),
                             y = c(s1 = 1, s2 = 0, s3 = 0, s4 = 1, s5 = 0, s6 = 1),
                             x = matrix(rep(1:6, 2), ncol = 2),
                             family = "binomial",
                             classify.rule = 0.5)

               )
})
