test_that("Pseudo-likelihood computing is correct", {
  expect_equal(pl_mrf2d(matrix(c(0,0,2,2), nrow = 2), mrfi(), theta_potts), log(exp(0)/(exp(0)+exp(0)+exp(-1)))*4, tolerance = 10^-6)
  expect_equal(pl_mrf2d(Z_potts, mrfi(), theta_potts*0),
               -log(3)*prod(dim(Z_potts)), tolerance = 10^-6)
  expect_is(pl_mrf2d(Z_potts, mrfi(), theta_potts*0, log_scale = FALSE), "numeric")

  # 'onepar'
  onepar <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "onepar")
  expect_setequal(names(onepar), c("theta", "value"))
  onepar <- fit_pl(Z_potts[1:30, 1:30], mrfi(), optim_args = list(method = "CG"), return_optim = TRUE)
  expect_true(all(c("theta", "value", "opt.convergence") %in% names(onepar)))
  expect_is(fit_pl(Z_potts[1:30,1:30], mrfi(), init = 0.3), "list")
  expect_error(fit_pl(Z_potts[1:30,1:30], mrfi(), init = c(0.3,0.3)))

  # 'oneeach'
  expect_is(fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "oneeach"), "list")

  # 'absdif'
  expect_is(fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "absdif"), "list")

  # 'dif'
  dif <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "dif")
  expect_setequal(names(dif), c("theta", "value"))

  # 'free'
  expect_is(fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "free"), "list")

  # with subregions
  Z2 <- ifelse( (row(Z_potts) + col(Z_potts)) %% 2, Z_potts, NA)
  expect_equal(pl_mrf2d(Z2, mrfi(), theta_potts), log(1/3)*length(Z2[!is.na(Z2)]))

  Z3 <- ifelse( row(Z_potts) < 120, Z_potts, NA)
  expect_is(fit_pl(Z3, mrfi(), "dif"), "list")
})
