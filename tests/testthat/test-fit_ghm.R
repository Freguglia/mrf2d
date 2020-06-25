test_that("fit_ghm works", {
  set.seed(1)
  Y <- matrix(rep(c(0,1,2), each = 300), 30, 30, byrow = TRUE) + rnorm(900, sd = 0.2)
  expect_is(fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 3) ,"hmrfout")
  expect_is(fit_ghm(Y, mrfi(), theta_potts, fixed_fn = polynomial_2d(c(1,0), dim(Y)), verbose = FALSE) ,"hmrfout")
  expect_is(fit_ghm(Y, mrfi(), theta_potts, init_mus = 0:2, init_sigmas = rep(0.2, 3), verbose = FALSE) ,"hmrfout")

  # Different Variances
  ## No fixed effect
  fit1 <- fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 3, equal_vars = FALSE)
  expect_true(all(abs(fit1$par$mu - 0:2) < 0.1))
  expect_true(all(abs(fit1$par$sigma - 0.2) < 0.1))

  ## With fixed effect
  fit2 <- fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 3, equal_vars = FALSE, fixed_fn = polynomial_2d(c(0,2), dim(Y)))
  expect_true(all(abs(fit2$par$mu - 0:2) < 0.1))
  expect_true(all(abs(fit2$par$sigma - 0.2) < 0.1))

  expect_equal(sd(fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 3, equal_vars = TRUE)$par$sigma),
               0)
  })

test_that("fit_ghm works on subregions", {
  set.seed(1)
  Y <- matrix(rep(c(0,1,2), each = 300), 30, 30, byrow = TRUE) + rnorm(900, sd = 0.2)
  Y <- ifelse(abs(row(Y) + col(Y)) <= 55, Y, NA)
  expect_is(fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 3) ,"hmrfout")
  expect_is(fit_ghm(Y, mrfi(), theta_potts, fixed_fn = polynomial_2d(c(2,0), dim(Y)), verbose = FALSE) ,"hmrfout")
  expect_is(fit_ghm(Y, mrfi(), theta_potts, init_mus = 0:2, init_sigmas = rep(0.2, 3), verbose = FALSE) ,"hmrfout")

  # Different variances
  ## No fixed effect
  fit1 <- fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 10, equal_vars = FALSE)
  expect_true(all(abs(fit1$par$mu - 0:2) < 0.1))
  expect_true(all(abs(fit1$par$sigma - 0.2) < 0.1))

  # With fixed effect
  fit2 <- fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 10, equal_vars = FALSE, fixed_fn = polynomial_2d(c(0,2), dim(Y)))
  expect_true(all(abs(fit2$par$mu - 0:2) < 0.1))
  expect_true(all(abs(fit2$par$sigma - 0.2) < 0.1))


  expect_equal(sd(fit_ghm(Y, mrfi(), theta_potts, verbose = FALSE, maxiter = 3, equal_vars = TRUE)$par$sigma),
               0)
})

