test_that("mrfi plotting works", {
  expect_is(mrfi(), "mrfi")
  expect_is(plot(mrfi()), "ggplot")
  expect_is(plot(mrfi(), include_axis = TRUE), "ggplot")
  expect_is(plot(mrfi(), include_opposite = FALSE), "ggplot")
})

test_that("mrfi creation works", {
  expect_error(mrfi(max_norm = -1))
  expect_error(mrfi(positions = c(1,0)))
  expect_error(mrfi(positions = list(c("1","0"))))
  expect_error(mrfi(positions = list( c(1.2, 1.5))))
  expect_identical(mrfi(1), mrfi(1, positions = list()))
  expect_identical(mrfi(1), rpositions(list(c(1,0), c(0,1))))
  expect_equal(mrfi(positions = list(c(3,3)))@Rmat, rbind(diag(2), c(3,3)))
})

test_that("mrfi validity checks work", {
  expect_true(mrfi_is_valid(mrfi(1)))
  m1 <- m2 <- mrfi(1)
  m1@Rmat <- matrix(0, 3, 3)
  expect_false(isTRUE(mrfi_is_valid(m1)))

  m2@Rmat <- m2@Rmat + 0.23
  expect_false(isTRUE(mrfi_is_valid(m2)))
})

test_that("mrfi subsetting and conversion", {
  expect_error(mrfi(1)[3])
  expect_error(mrfi(1)[[3]])
  expect_identical(mrfi(1), mrfi(1)[1:2])
  expect_identical(as.list(mrfi(1)), mrfi(1)[[1:2]])

  expect_equal(length(mrfi(1)), 2)

  expect_identical(mrfi_to_string(mrfi(1)), "{(1,0),(0,1)}")
})

test_that("mrfi + and - operators work", {
  expect_error(mrfi(1) + c(1,0,0))
  expect_error(mrfi(1) + 2)
  expect_error(mrfi(1) + c(1.2, 0))
  expect_identical(mrfi(1, positions = list(c(2,0))), mrfi(1) + c(2,0))
  expect_identical(mrfi(1, positions = list(c(2,0))), mrfi(1) + mrfi(0, positions = list(c(2,0))))
  expect_identical(mrfi(1), mrfi(1) + mrfi(1))
  expect_identical(mrfi(1), mrfi(1) + c(-1,0))

  expect_error(mrfi(1) - c(1,0,0))
  expect_error(mrfi(1) - 2)
  expect_error(mrfi(1) - c(1.2,0))
  expect_identical(mrfi(1), mrfi(1) - c(2,0))
  expect_identical(mrfi(1), mrfi(1) - c(0,0))
  expect_identical(mrfi(0, positions = list(c(0,1))), mrfi(1) - c(1,0))
  expect_identical(mrfi(1) - c(1,0), mrfi(1) - c(-1,0))
  expect_identical(mrfi(1) - mrfi(1), mrfi(0))
  expect_setequal(as.list(mrfi(2) - mrfi(1)), as.list(mrfi(0, positions = list(c(2,0), c(0,2), c(1,1), c(-1,1)))))
  expect_identical(mrfi(1) - mrfi(0, positions = list(c(1,0))), mrfi(1) - mrfi(0, positions = list(c(-1,0))))
})
