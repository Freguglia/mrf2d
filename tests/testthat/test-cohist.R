test_that("cohist works", {
  expect_identical(dim(cohist(Z_potts, mrfi(1))), c(3L,3L,2L))
})
