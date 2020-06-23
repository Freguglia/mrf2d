test_that("mrfout prints correctly", {
  # Print
  a <- fit_pl(Z_potts[1:30, 1:30], mrfi(1))
  prnt <- capture.output(print(a))
  expect_true(sum(grepl("Pseudolikelihood", prnt)) > 0)

  # Summary
  a <- fit_pl(Z_potts[1:30, 1:30], mrfi(1), family = "oneeach")
  prnt <- capture.output(summary(a))
  expect_true(sum(grepl("Interaction", prnt)) > 0)

  a <- fit_pl(Z_potts[1:30, 1:30], mrfi(1), family = "absdif")
  prnt <- capture.output(summary(a))
  expect_true(sum(grepl("Interaction", prnt)) > 0)

  a <- fit_pl(Z_potts[1:30, 1:30], mrfi(1), family = "dif")
  prnt <- capture.output(summary(a))
  expect_true(sum(grepl("Interaction", prnt)) > 0)

  a <- fit_pl(Z_potts[1:30, 1:30], mrfi(1), family = "free")
  prnt <- capture.output(summary(a))
  expect_true(sum(grepl("Interaction", prnt)) > 0)
})
