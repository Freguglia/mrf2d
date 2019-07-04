test_that("mrfi class works", {
  expect_is(mrfi(), "mrfi")
  expect_is(plot(mrfi()), "ggplot")
  expect_is(plot(mrfi(), no_axis = TRUE), "ggplot")

  ## Creator tests
  expect_error(mrfi(max_norm = -1))
  expect_error(mrfi(positions = c(1,0)))
  expect_error(mrfi(positions = list(c("1","0"))))
  expect_identical(mrfi(positions = list(c(3,3)))@Rmat, rbind(diag(2), c(3,3)))
})
