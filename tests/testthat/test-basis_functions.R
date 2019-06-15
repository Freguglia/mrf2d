test_that("fourier_2d", {
  expect_error(fourier_2d(1, c(3,4)))
  expect_error(fourier_2d(c(-1,4), c(5,5)))
  expect_error(fourier_2d(c(1,4), c(15,15,2)))
})
