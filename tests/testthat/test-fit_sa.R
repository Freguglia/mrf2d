test_that("fit_sa works", {
  expect_error(fit_sa(Z_potts, mrfi(), family = "onepar", init = c(1,2), gamma_seq = 1:0))
  expect_is(fit_sa(Z_potts, mrfi(), family = "onepar", gamma_seq = 1:0), "list")
  expect_true(is_valid_array(fit_sa(Z_potts, mrfi(), family = "dif", init = rep(-1, 4*2), gamma_seq = seq(1,0,-0.1))$theta, "dif"))
  expect_is(fit_sa(Z_potts, mrfi(), family = "onepar", gamma_seq = seq(1,0,-0.1), refresh_each = 2, refresh_cycles = 5), "list")
  })

test_that("fit_sa works with subregions", {
  Z <- Z_potts
  Z <- ifelse( col(Z) >= (row(Z)- 75)^2/150, Z, NA )
  expect_error(fit_sa(Z, mrfi(), family = "onepar", init = c(1,2), gamma_seq = 1:0))
  expect_is(fit_sa(Z, mrfi(), family = "onepar", gamma_seq = 1:0), "list")
  expect_true(is_valid_array(fit_sa(Z, mrfi(), family = "dif", init = rep(-1, 4*2), gamma_seq = seq(1,0,-0.1))$theta, "dif"))
  expect_is(fit_sa(Z, mrfi(), family = "onepar", gamma_seq = seq(1,0,-0.1), refresh_each = 2, refresh_cycles = 5), "list")
})
