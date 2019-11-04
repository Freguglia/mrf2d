test_that("mrfi class works", {
  expect_is(mrfi(), "mrfi")
  expect_is(plot(mrfi()), "ggplot")
  expect_is(plot(mrfi(), no_axis = TRUE), "ggplot")

  ## Creator tests
  expect_error(mrfi(max_norm = -1))
  expect_error(mrfi(positions = c(1,0)))
  expect_error(mrfi(positions = list(c("1","0"))))
  expect_identical(mrfi(positions = list(c(3,3)))@Rmat, rbind(diag(2), c(3,3)))

  ## Subsetting and conversion
  expect_error(mrfi(1)[3])
  expect_error(mrfi(1)[[3]])
  expect_identical(mrfi(1), mrfi(1)[1:2])
  expect_identical(as.list(mrfi(1)), mrfi(1)[[1:2]])

  # '+'
  expect_error(mrfi(1) + c(1,0,0))
  expect_error(mrfi(1) + 2)
  expect_identical(mrfi(1, positions = list(c(2,0))), mrfi(1) + c(2,0))
  expect_identical(mrfi(1, positions = list(c(2,0))), mrfi(1) + mrfi(0, positions = list(c(2,0))))
})
