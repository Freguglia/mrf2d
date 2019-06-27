test_that("mrfi class works", {
  expect_is(ipotts, "mrfi")
  expect_is(plot(ipotts), "ggplot")
  expect_warning(ipotts, regexp = NA)
})
