test_that("hmrfout methods work", {
  m <- mrfi(1) + c(4,4)
  basis <- polynomial_2d(c(1,1), dim(hfield1))
  a <- fit_ghm(hfield1, m, expand_array(c(-1,-1,.2), "oneeach", m, 1),
               fixed_fn = basis, init_mus = c(5,9), init_sigmas = c(1.3,1.3),
               maxiter = 2)
  y <- hfield1
  y[3,] <- NA
  b <- fit_ghm(hfield1, m, expand_array(c(-1,-1,.2), "oneeach", m, 1),
               fixed_fn = basis, init_mus = c(5,9), init_sigmas = c(1.3,1.3),
               maxiter = 2)

  prnta <- capture.output(print(a))
  prntb <- capture.output(print(b))
  expect_true(sum(grepl("EM", prnta)) > 0)
  expect_true(sum(grepl("EM", prntb)) > 0)

  smra <- capture.output(summary(a))
  smrb <- capture.output(summary(b))
  expect_true(sum(grepl("mu", smra)) > 0)
  expect_true(sum(grepl("mu", smrb)) > 0)

  pdf(NULL)
  expect_is(plot(a), "hmrfout")
  expect_is(plot(b), "hmrfout")
})
