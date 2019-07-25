test_that("Gibbs Sampler works", {
  expect_error(rmrf2d("string", mrfi(), theta_potts))
  expect_error(rmrf2d(c(20,20,20), mrfi(), theta_potts))
  expect_error(rmrf2d(c(-30,30), mrfi(), theta_potts))

  set.seed(1)
  Z <- rmrf2d(c(30,30), mrfi(), theta_potts)

  expect_type(Z, "integer")
  expect_true(is.matrix(Z))
  expect_true(all(Z %in% 0:3))
})
