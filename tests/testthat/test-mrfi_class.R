test_that("mrfi class works", {
  expect_is(ipotts, "mrfi")
  expect_is(plot(ipotts), "ggplot")
  expect_is(plot(ipotts, no_axis = TRUE), "ggplot")
})
