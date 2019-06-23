test_that("mrf Gibbs Sampler", {
  theta <- array(diag(3)*2 -1, dim = c(3,3,2))*0.7
  potts <- new("mrfi", Rmat = diag(2), n_neis = 2)

  expect_error(rmrf2d("string", potts, theta))
  expect_error(rmrf2d(c(20,20,20), potts, theta))
  expect_error(rmrf2d(c(-30,30), potts, theta))

  set.seed(1)
  Z <- rmrf2d(c(30,30), potts, theta)

  expect_type(Z, "integer")
  expect_true(is.matrix(Z))
  expect_true(all(Z %in% 0:3))
})
