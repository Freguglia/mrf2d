test_that("Conditional probabilities are correct", {
  z <- matrix(c(1,0,0,
                0,1,0,
                0,0,0), nrow = 3, byrow = TRUE)
  theta <- -array(1 - diag(2), dim = c(2,2,2))

  expect_equivalent(cp_mrf2d(z, mrfi(1), theta, pos = c(1,1)),
               c(1/(1+exp(-2)), exp(-2)/(1 + exp(-2))))
  expect_equivalent(cp_mrf2d(z, mrfi(1), theta, pos = c(2,2)),
                    c(1/(1+exp(-4)), exp(-4)/(1 + exp(-4))))

  zi <- z
  zi[3,2] <- NA
  expect_equivalent(cp_mrf2d(zi, mrfi(1), theta, pos = c(1,1)),
                    c(1/(1+exp(-2)), exp(-2)/(1 + exp(-2))))
  expect_equivalent(cp_mrf2d(zi, mrfi(1), theta, pos = c(2,2)),
                    c(1/(1+exp(-3)), exp(-3)/(1 + exp(-3))))

})
