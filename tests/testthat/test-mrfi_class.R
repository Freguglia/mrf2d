test_that("mrfi class works", {
  expect_is(mrfi(), "mrfi")
  expect_is(plot(mrfi()), "ggplot")
  expect_is(plot(mrfi(), no_axis = TRUE), "ggplot")
})
