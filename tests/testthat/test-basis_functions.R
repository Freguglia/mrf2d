test_that("fourier_2d", {
  expect_error(fourier_2d(max_freqs = 1, lattice_dim = c(3,4)))
  expect_error(fourier_2d(max_freqs = c(-1,4), lattice_dim = c(5,5)))
  expect_error(fourier_2d(max_freqs = c(1,4), lattice_dim = c(15,15,2)))
  expect_is(fourier_2d(max_freqs = c(2,2), lattice_dim = c(10,10)), "list")
  expect_is(fourier_2d(max_freqs = c(2,2), lattice_dim = c(10,10))[[1]], "function")
})

test_that("polynomial_2d", {
  expect_error(polynomial_2d(poly_deg = 1, lattice_dim = c(3,4)))
  expect_error(polynomial_2d(poly_deg = c(-1,4), lattice_dim = c(5,5)))
  expect_error(polynomial_2d(poly_deg = c(1,4), lattice_dim = c(15,15,2)))
  expect_is(polynomial_2d(poly_deg = c(2,2), lattice_dim = c(10,10)), "list")
  expect_is(polynomial_2d(poly_deg = c(2,2), lattice_dim = c(10,10))[[1]], "function")
})
