test_that("fit_icm works", {
  expect_is(fit_icm(Z_potts, mrfi(), theta_potts, error_prob = 0.2, 6), "matrix")
  expect_error(fit_icm(Z_potts, mrfi(), theta_potts, error_prob = 0, 6))
})
