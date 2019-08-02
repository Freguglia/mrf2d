test_that("fit_sa works", {
  expect_error(fit_sa(Z_potts, mrfi(), family = "onepar", init = c(1,2), gamma_seq = 1:0))
  expect_is(fit_sa(Z_potts, mrfi(), family = "onepar", gamma_seq = 1:0), "list")
  expect_true(is_valid_array(fit_sa(Z_potts, mrfi(), family = "dif", init = rep(-1, 4*2), gamma_seq = seq(1,0,-0.05))$theta, "dif"))
})