test_that("Pseudo-likelihood computing is correct", {
  expect_equal(pl_mrf2d(matrix(c(0,0,2,2), nrow = 2), ipotts, theta_potts), log(exp(0)/(exp(0)+exp(0)+exp(-1)))*4, tolerance = 10^-6)
  expect_equal(pl_mrf2d(Z_potts, ipotts, theta_potts*0),
               -log(3)*prod(dim(Z_potts)), tolerance = 10^-6)

  # 'onepar' maximum pseudo-likelihood fitting
  onepar <- fit_pl(Z_potts[1:30, 1:30], ipotts, family = "onepar")
  expect_setequal(names(onepar), c("theta", "value"))
  onepar <- fit_pl(Z_potts[1:30, 1:30], ipotts, optim_args = list(method = "CG"), return_optim = TRUE)
  expect_true(all(c("theta", "value", "opt.convergence") %in% names(onepar)))

  # 'dif' maximum pseudo-likelihood fitting
  dif <- fit_pl(Z_potts[1:30, 1:30], ipotts, family = "dif")
  expect_setequal(names(dif), c("theta", "value"))
})
