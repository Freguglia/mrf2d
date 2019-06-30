test_that("vector - array conversions works", {
  # 'onepar' family
  expect_identical(vec_to_array(-1, "onepar", C = 2, n_R = 2), theta_potts)
  expect_identical(array_to_vec(vec_to_array(-1, "onepar", C = 2, n_R = 2), "onepar"), -1)
  expect_identical(array_to_vec(vec_to_array(-1, "onepar", C = 6, n_R = 4), "onepar"), -1)

  # 'oneeach' family
  expect_identical(array_to_vec(vec_to_array(c(1,2), "oneeach", C = 2, n_R = 2), "oneeach"), c(1,2))
  expect_identical(array_to_vec(vec_to_array(c(1,2,3,4), "oneeach", C = 6, n_R = 4), "oneeach"), c(1,2,3,4))
  expect_error(vec_to_array(c(1,2), "oneeach", C = 2, n_R = 3))
  expect_error(vec_to_array(c(1,2), "onepar", C = 2, n_R = 2))

  # 'absdif' family
  set.seed(1); v1 <- runif(2*2); v2 <- runif(3*3)
  expect_identical(array_to_vec(vec_to_array(v1, "absdif", C = 2, n_R = 2), "absdif"), v1)
  expect_identical(array_to_vec(vec_to_array(v2, "absdif", C = 3, n_R = 3), "absdif"), v2)
  expect_equal(vec_to_array(v1, "absdif", 2, 2)[3,1,1], v1[2])

  # 'dif' family
  set.seed(1); v1 <- runif(4*2); v2 <- runif(6*3)
  expect_identical(array_to_vec(vec_to_array(v1, "dif", C = 2, n_R = 2), "dif"), v1)
  expect_identical(array_to_vec(vec_to_array(v2, "dif", C = 3, n_R = 3), "dif"), v2)
  expect_equal(vec_to_array(v1, "dif", 2, 2)[3,1,1], v1[1])

})
