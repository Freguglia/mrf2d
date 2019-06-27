test_that("plots work", {
  expect_is(dplot(Z_potts), "ggplot")
  expect_is(cplot(Z_potts + rnorm(prod(dim(Z_potts)))), "ggplot")
})
