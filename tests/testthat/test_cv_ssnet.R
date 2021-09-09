test_that("Bad inputs are not accepted.", {
  expect_error(cv_ssnet(model = "ss",
                        x = matrix(rnorm(10*10), nrow = 10, ncol = 10),
                        y = rnorm(10),
                        foldid = list(f1 = sample(1:3, 10, T),
                                      f2 = sample(1:3, 10, T))),
               "If foldid is supplied, must be vector, matrix, or data frame."
               )
  expect_error(cv_ssnet(model = "ss",
                        x = matrix(rnorm(10*10), nrow = 10, ncol = 10),
                        y = rnorm(10),
                        foldid = data.frame(
                          f1 = sample(1:3, 10, T),
                          f2 = sample(1:4, 10, T))),
               "Each column of foldid must have the same number of unique folds.")

})
