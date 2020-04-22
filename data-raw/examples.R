# Example binary field with sparse interaction structure.
library(mrf2d)
set.seed(1)
m1 <- mrfi(1) + c(4,4)
theta1 <- mrf2d:::vec_to_array(c(-1, -1, .2), family = "oneeach", C = 1, 3)
field1 <- rmrf2d(c(150,150), m1, theta1, cycles = 80)

hfield1 <- (field1 + 3) + 
  rnorm(prod(dim(field1)), mean = 0, sd = (as.vector(field1) + 1)*0.6) +
  0.05*col(field1)

usethis::use_data(field1, hfield1, overwrite = TRUE)
