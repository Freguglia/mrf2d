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

test_that("Gibbs Sampler works with sub_region and fixed_region", {
  set.seed(1)
  mask <- matrix(TRUE, 30, 30)
  mask <- ifelse((row(mask)-15)^2 + (col(mask)-15)^2 <= 15^2, TRUE, FALSE)

  masked <- rmrf2d(c(30,30), mrfi(), theta_potts, 5, sub_region = mask)
  expect_error(rmrf2d(c(30,30), mrfi(), theta_potts, 5, sub_region = c(30,30)))

  expect_is(masked, "matrix")
  expect_is(rmrf2d(Z_potts[1:30, 1:30], mrfi(), theta_potts, 5, sub_region = mask), "matrix")
  expect_true(any(is.na(rmrf2d(masked, mrfi(), theta_potts, 5))))

  expect_true(any(is.na(masked)))

  expect_error(rmrf2d(c(35,35), mrfi(), theta_potts, 5, sub_region = mask))
  expect_error(rmrf2d(c(25,25), mrfi(), theta_potts, 5, sub_region = mask))
  expect_error(rmrf2d(masked_noNA[1:25,1:25], mrfi(), theta_potts, 5, sub_region = mask))
  expect_warning(rmrf2d(masked, mrfi(), theta_potts, 5, sub_region = mask), "'init_Z' has NA values")

  # Only fixed
  border <- matrix(FALSE, 30, 31)
  border[c(1,30),] <- border[,c(1,31)] <- TRUE
  init <- matrix(sample(0:2, replace = TRUE, 30*31), 30, 31)
  init[border] <- 0
  expect_error(rmrf2d(init, mrfi(), theta_potts, 5, fixed_region = border[1:30,1:30]))
  expect_is(rmrf2d(init, mrfi(), theta_potts, 5, fixed_region = border*1), "matrix")
  expect_error(rmrf2d(c(30,30), mrfi(), theta_potts, 5, fixed_region = border))
  expect_is(rmrf2d(c(30,31), mrfi(), theta_potts, 40, fixed_region = border), "matrix")

  # Both fixed and sub-regions
  init_diag <- sample(0:2, 30*30, replace = TRUE)
  init_diag <- init_diag*(1-ifelse(diag(30), TRUE, FALSE))
  expect_warning(rmrf2d(init_diag, mrfi(), theta_potts, 40, fixed_region = ifelse(diag(30),TRUE,FALSE), sub_region = mask), "Some pixels")
  expect_warning(rmrf2d(c(30,30), mrfi(), theta_potts, 5, fixed_region = ifelse(diag(30),TRUE,FALSE), sub_region = mask), "Some pixels")
})
