test_that("smr_array works", {
  expect_length(smr_array(theta_potts, "onepar"), 1)
  expect_length(smr_array(theta_potts, "oneeach"), 2)
  expect_length(smr_array(theta_potts, "absdif"), 2*2)
  expect_length(smr_array(theta_potts, "dif"), 2*2*2)
  expect_length(smr_array(theta_potts, "free"), 2*(3^2 - 1))
})

test_that("smr_stat works", {
  expect_length(smr_stat(Z_potts, mrfi(), "onepar"), 1)
  expect_length(smr_stat(Z_potts, mrfi(), "oneeach"), 2)
  expect_length(smr_stat(Z_potts, mrfi(), "absdif"), 2*2)
  expect_length(smr_stat(Z_potts, mrfi(), "dif"), 2*2*2)
  expect_length(smr_stat(Z_potts, mrfi(), "free"), 2*(3^2 - 1))
})

