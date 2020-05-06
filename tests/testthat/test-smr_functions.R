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

test_that("expand_array", {
  expect_is(expand_array(0.95, family = "onepar", mrfi(1), C = 1), "array")
  expect_is(expand_array(c(0.1, 0.2), family = "oneeach", mrfi(1), C = 1), "array")

  expect_equivalent(1:2, smr_array(expand_array(1:2, "oneeach", mrfi(1), 2), "oneeach"))
  expect_equivalent(1:4, smr_array(expand_array(1:4, "absdif", mrfi(1), 2), "absdif"))
  expect_equivalent(1:8, smr_array(expand_array(1:8, "dif", mrfi(1), 2), "dif"))
  expect_equivalent(1:16, smr_array(expand_array(1:16, "free", mrfi(1), 2), "free"))
})

