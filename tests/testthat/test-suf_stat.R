test_that("sufficient statistics computing works", {
  pc <- table_relative_3d(Z_potts, diag(2), C = 2)
  expect_equal(length(suf_stat(pc, "onepar")), 1)
  expect_equal(length(suf_stat(pc, "oneeach")), 2)
  expect_equal(length(suf_stat(pc, "absdif")), 2*2)
  expect_equal(length(suf_stat(pc, "dif")), 2*4)
  expect_equal(length(suf_stat(pc, "free")), 2*(3*3 -1))
})
