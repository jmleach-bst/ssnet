test_that("Check formating for IAR priors", {
  expect_error(ssnet(x = matrix(rnorm(100), nrow = 10, ncol = 10),
                     y = rnorm(10), family = "gaussian",
                     iar.prior = TRUE),
               "User must specify both im.res and x.")
  expect_error(ssnet(x = matrix(rnorm(10 * 16), nrow = 10, ncol = 16),
                     y = rnorm(10), family = "gaussian",
                     iar.prior = TRUE, im.res = c(5, 5)),
               "The product of im.res should equal the number of columns in x.")
  expect_error(format_iar(x = matrix(rnorm(100), nrow = 10, ncol = 10)),
               "User must specify both im.res and x.")
  expect_error(format_iar(x = matrix(rnorm(10 * 16), nrow = 10, ncol = 16),
                          im.res = c(5, 5)),
               "The product of im.res should equal the number of columns in x.")
  expect_error(format_iar(x = matrix(rnorm(10 * 16), nrow = 10, ncol = 16),
                          im.res = c(4, 4), tau.prior = "manual"),
               "User must specify tau.manual")
  expect_type(format_iar(x = matrix(rnorm(10 * 16), nrow = 10, ncol = 16),
                         im.res = c(4, 4)),
              "list")

})

test_that("Check Q2 maximization output", {
  expect_type(max_q2_iar(iar.data = format_iar(x = matrix(rnorm(10 * 16),
                                                          nrow = 10, ncol = 16),
                                               im.res = c(4, 4)),
                         p = rnorm(16)), "double")
  expect_equal(length(max_q2_iar(iar.data = format_iar(x = matrix(rnorm(10 * 16),
                                                                  nrow = 10, ncol = 16),
                                                       im.res = c(4, 4)),
                                 p = rnorm(16))),
               length(rnorm(16)))
})
