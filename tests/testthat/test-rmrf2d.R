test_that("mrf Gibbs Sampler", {
  expect_error(rmrf2d("string", ipotts, theta_potts))
  expect_error(rmrf2d(c(20,20,20), ipotts, theta_potts))
  expect_error(rmrf2d(c(-30,30), ipotts, theta_potts))

  set.seed(1)
  Z <- rmrf2d(c(30,30), ipotts, theta_potts)

  expect_type(Z, "integer")
  expect_true(is.matrix(Z))
  expect_true(all(Z %in% 0:3))
})
