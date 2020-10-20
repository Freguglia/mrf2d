test_that("rmrf2d_mc works as expected", {
  set.seed(1)
  a <- rmrf2d_mc(c(80, 80), mrfi(1), theta_potts, family = "oneeach", burnin = 3,
                 nmc = 10)
  expect_equal(nrow(a), 10)
  expect_is(a, "matrix")
})
