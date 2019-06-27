test_that("mrfi class works", {
  expect_is(ipotts, "mrfi")
  expect_is(plot(ipotts), "ggplot")
  expect_is(plot(ipotts, no_axis = TRUE), "ggplot")
  expect_warning(show(ipotts), regexp = NA)
  expect_warning(print(ipotts), regexp = NA)
})
