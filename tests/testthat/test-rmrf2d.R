test_that("Gibbs Sampler works", {
  expect_error(rmrf2d("string", mrfi(), theta_potts))
  expect_error(rmrf2d(c(20,20,20), mrfi(), theta_potts))
  expect_error(rmrf2d(c(-30,30), mrfi(), theta_potts))

  set.seed(1)
  Z <- rmrf2d(c(30,30), mrfi(), theta_potts)

  expect_type(Z, "integer")
  expect_true(is.matrix(Z))
  expect_true(all(Z %in% 0:3))

  t2 <- array(theta_potts, dim = c(3,3,3))
  t2[,,3] <- 0

  expect_is(rmrf2d(c(30,30), mrfi(1, positions = list(c(3,3))),t2), "matrix")

})

test_that("Gibbs Sampler works with sub-lattices", {
  set.seed(1)
  mask <- matrix(TRUE, 30, 30)
  mask <- ifelse((row(mask)-15)^2 + (col(mask)-15)^2 <= 15^2, TRUE, FALSE)

  masked <- rmrf2d(c(30,30), mrfi(), theta_potts, 40, sub_lattice = mask)
  expect_error(rmrf2d(c(30,30), mrfi(), theta_potts, 40, sub_lattice = c(30,30)))

  expect_is(masked, "matrix")
  expect_is(rmrf2d(Z_potts[1:30, 1:30], mrfi(), theta_potts, 40, sub_lattice = mask), "matrix")
  expect_true(any(is.na(rmrf2d(masked, mrfi(), theta_potts, 40))))

  masked_noNA <- rmrf2d(c(30,30), mrfi(), theta_potts, 40, sub_lattice = mask, mask_na = FALSE)
  expect_is(masked_noNA, "matrix")

  expect_true(any(is.na(masked)))
  expect_false(any(is.na(masked_noNA)))

  expect_error(rmrf2d(c(35,35), mrfi(), theta_potts, 40, sub_lattice = mask))
  expect_error(rmrf2d(c(25,25), mrfi(), theta_potts, 40, sub_lattice = mask))
  expect_error(rmrf2d(masked_noNA[1:25,1:25], mrfi(), theta_potts, 40, sub_lattice = mask))
  expect_warning(rmrf2d(masked, mrfi(), theta_potts, 40, sub_lattice = mask), "'init_Z' has NA values")
})
