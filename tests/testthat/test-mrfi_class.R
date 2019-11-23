test_that("mrfi plotting works", {
  expect_is(mrfi(), "mrfi")
  expect_is(plot(mrfi()), "ggplot")
  expect_is(plot(mrfi(), no_axis = TRUE), "ggplot")
})

test_that("mrfi creation works", {
  expect_error(mrfi(max_norm = -1))
  expect_error(mrfi(positions = c(1,0)))
  expect_error(mrfi(positions = list(c("1","0"))))
  expect_identical(mrfi(1), mrfi(1, positions = list()))
  expect_equal(mrfi(positions = list(c(3,3)))@Rmat, rbind(diag(2), c(3,3)))
})

test_that("mrfi subsetting and conversion", {
  expect_error(mrfi(1)[3])
  expect_error(mrfi(1)[[3]])
  expect_identical(mrfi(1), mrfi(1)[1:2])
  expect_identical(as.list(mrfi(1)), mrfi(1)[[1:2]])
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
