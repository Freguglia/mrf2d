test_that("Pseudo-likelihood computing is correct", {
  theta <- array(diag(3)*2 -1, dim = c(3,3,2))*0.7
  potts <- new("mrfi", Rmat = diag(2), n_neis = 2)
  set.seed(1)
  Z <- rmrf2d(c(30,30), potts, theta)
  expect_equal(pl_mrf2d(matrix(c(0,0,2,2), nrow = 2), potts, theta), log(exp(0)/(exp(0)+exp(0)+exp(-1.4)))*4, tolerance = 10^-6)
  expect_equal(pl_mrf2d(Z, potts, theta*0), -log(3)*30*30, tolerance = 10^-6)
})
