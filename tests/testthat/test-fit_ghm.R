test_that("fit_ghm works", {
  set.seed(1)
  Y <- Z_potts + rnorm(length(Z_potts), sd = 0.2)
  expect_is(fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE) ,"list")
  expect_is(fit_ghm(Y, mrfi(), theta_potts, fixed_fn = polynomial_2d(c(2,2), dim(Y)), verbose = FALSE) ,"list")
  expect_is(fit_ghm(Y, mrfi(), theta_potts, init_mus = 0:2, init_sigmas = rep(0.2, 3), verbose = FALSE) ,"list")
})
