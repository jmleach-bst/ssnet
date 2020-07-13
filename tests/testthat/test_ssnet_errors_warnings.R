test_that("Correct errors and warnings are raised for bad input", {
  # bad x and/or y inputs
  expect_error(ssnet(x = list(a = 'x'), family = "gaussian"), "x should be a matrix")
  expect_error(ssnet(x = data.frame(a = 'x'), family = "gaussian"), "x should be a matrix")
  expect_error(ssnet(x = c(a = 'x'), family = "gaussian"), "x should be a matrix")
  expect_error(ssnet(x = matrix(1:100, nrow = 10, ncol = 10), y = 1,
                     family = "gaussian"), "length of y should equal number of rows in x")
  expect_error(ssnet(x = matrix(1:100, nrow = 10, ncol = 10), y = rep("a", 10),
                         family = "gaussian"),
               "gaussian family requires numeric y")
  expect_error(ssnet(x = matrix(1:100, nrow = 10, ncol = 10), y = rnorm(10),
                     init = rep(0, 5),
                     family = "gaussian"),
               "must specify initial value to each coefficient except intercept")

  # bad EN parameter inputs
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                         y = rnorm(10), family = "gaussian",
                         alpha = -1),
               "alpha must be 0 or greater and cannot exceed 1.")
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                         y = rnorm(10), family = "gaussian",
                         alpha = 1.01),
               "alpha must be 0 or greater and cannot exceed 1.")

  # bad spike/slab scales
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                     y = rnorm(10), family = "gaussian",
                     ss = c(-1, 1)),
               "scale values must exceed 0.")
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                     y = rnorm(10), family = "gaussian",
                     ss = c(0, 1)),
               "scale values must exceed 0.")
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                     y = rnorm(10), family = "gaussian",
                     ss = c(0.05, 1, 2)),
               "ss should contain only 2 elements, 1 spike and 1 slab scale.")
})
