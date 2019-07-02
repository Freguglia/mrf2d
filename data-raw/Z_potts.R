## code to prepare `Z_potts` dataset goes here
set.seed(1)
theta_potts <- mrf2d:::vec_to_array(-1, "onepar", 2, 2)
Z_potts <- rmrf2d(c(150,150), mrfi(), theta_potts, steps = 80)

usethis::use_data(Z_potts, mrfi(), theta_potts, overwrite = TRUE)
