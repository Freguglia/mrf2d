## code to prepare `Z_potts` dataset goes here
set.seed(1)
ipotts <- new("mrfi", Rmat = diag(2), n_neis = 2)
theta_potts <- mrf2d:::vec_to_array(0.7, "onepar", 2, 2)
Z_potts <- rmrf2d(c(150,150), potts_interaction, potts_theta, steps = 80)

usethis::use_data(Z_potts, ipotts, theta_potts, overwrite = TRUE)
