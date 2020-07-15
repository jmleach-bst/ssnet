test_that("Bad inputs are not accepted.", {
  expect_error(compare_ssnet(variable_selection = TRUE,
                             x = matrix(rnorm(100), nrow = 10, ncol = 10),
                             B = c(1, 2)),
               "Variable selection measures requires a true parameter vector of appropriate length. \n
            Did you forget to remove the intercept or specify B?")
  expect_error(compare_ssnet(models = "read my mind!"),
               "Models must be a combination of all, glmnet, ss, ss_iar.")
  expect_error(compare_ssnet(models = "glmnet", family = "read my mind!"),
               "Models must be a combination of gaussian, binomial, poisson, cox.")

})
