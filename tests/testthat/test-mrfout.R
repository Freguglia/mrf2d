test_that("mrfout prints correctly", {
  a <- fit_pl(Z_potts, mrfi(1))
  prnt <- capture.output(print(a))
  expect_true(sum(grepl("Pseudolikelihood", prnt)) > 0)
})
