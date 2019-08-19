# @name fit_icm
# @author Victor Freguglia
# @title Iterated Conditional Modes algorithm for image restoration
#
# @description Computes maximum probability configuration for a Markov random
#  with corrupted pixel values.
#
# @inheritParams pl_mrf2d
# @param n_cycles The number of times each pixel is updated.
# @param error_prob The probability a pixel's value is corrupted.
#
# @return A `matrix` with the resulting configuration.
#' @keywords internal
fit_icm <- function(init_Z, mrfi, theta, error_prob, cycles = 10){
  if(error_prob >= 1 | error_prob <= 0) {
    stop("Error probability 'error_prob' must be between 0 and 1.")}
  return(icm_restoration_cpp(init_Z, mrfi@Rmat, theta, error_prob, cycles))
}
