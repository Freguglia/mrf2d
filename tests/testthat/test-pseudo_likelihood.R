test_that("Pseudo-likelihood computing is correct", {
  expect_equal(pl_mrf2d(matrix(c(0,0,2,2), nrow = 2), mrfi(), theta_potts), log(exp(0)/(exp(0)+exp(0)+exp(-1)))*4, tolerance = 10^-6)
  expect_equal(pl_mrf2d(Z_potts, mrfi(), theta_potts*0),
               -log(3)*prod(dim(Z_potts)), tolerance = 10^-6)
  expect_is(pl_mrf2d(Z_potts, mrfi(), theta_potts*0, log_scale = FALSE), "numeric")
})

test_that("fit_pl is returning correctly", {
  # 'onepar'
  onepar <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "onepar")
  expect_setequal(names(onepar), c("theta", "mrfi", "family", "method", "value", "Z"))
  onepar <- fit_pl(Z_potts[1:30, 1:30], mrfi(), return_optim = TRUE)
  expect_true(all(c("theta", "value", "opt.convergence") %in% names(onepar)))
  expect_is(fit_pl(Z_potts[1:30,1:30], mrfi(), init = -0.8), "mrfout")
  expect_error(fit_pl(Z_potts[1:30,1:30], mrfi(), init = c(0.3,0.3)))
})

test_that("fit_pl works with family 'oneeach'", {
  oneeach <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "oneeach")
  expect_is(oneeach, "mrfout")
})

test_that("fit_pl works with family 'absdif'", {
  absdif <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "absdif")
  expect_is(absdif, "mrfout")
})

test_that("fit_pl works with family 'dif'", {
  dif <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "dif")
  expect_setequal(names(dif), c("theta", "mrfi", "family", "method", "value", "Z"))
})

test_that("fit_pl works with family 'free'", {
  free1 <- fit_pl(Z_potts[1:30, 1:30], mrfi(), family = "free")
  expect_is(free1, "mrfout")
})

test_that("Pseudo-likelihood works with sub_region", {
  Z2 <- ifelse( (row(Z_potts) + col(Z_potts)) %% 2, Z_potts, NA)
  expect_equal(pl_mrf2d(Z2, mrfi(), theta_potts), log(1/3)*length(Z2[!is.na(Z2)]))

  Z3 <- ifelse( (row(Z_potts) < 31) & (col(Z_potts) < 31), Z_potts, NA)
  expect_is(fit_pl(Z3, mrfi(), "dif"), "mrfout")
})
