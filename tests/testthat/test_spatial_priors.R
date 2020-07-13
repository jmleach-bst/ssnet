devtools::load_all()

test_that("Spatial priors are correctly applied", {
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                     y = rnorm(10), family = "gaussian",
                     iar.prior = TRUE),
               "Require image dimensions, im.res, to create adjacency matrix.")
  expect_error(ssnet(x = matrix(rnorm(10 * 16), nrow = 10, ncol = 16),
                     y = rnorm(10), family = "gaussian",
                     iar.prior = TRUE, im.res = c(5, 5)),
               "The product of im.res should equal the number of columns in x.")
  expect_error(format_iar(x = matrix(rnorm(100), nrow = 10, ncol = 10)),
               "Require image dimensions, im.res, to create adjacency matrix.")
  expect_error(format_iar(x = matrix(rnorm(10 * 16), nrow = 10, ncol = 16),
                          im.res = c(5, 5)),
               "The product of im.res should equal the number of columns in x.")
})
