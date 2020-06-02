#' @name cpmrf2d
#' @author Victor Freguglia
#' @title Conditional probabilities in a pixel position
#'
#' @description Computes the vector of conditional probabilities
#' for a pixel position given a field, an interaction structure and
#' a parameter array.
#'
#' @inheritParams pl_mrf2d
#' @param pos Length-2 vector with the position to compute conditional
#' probabilities in.
#'
#' @return A `numeric` vector with the conditional probabilities.
#'
#' @examples
#' cp_mrf2d(Z_potts, mrfi(1), theta_potts, c(57,31))
#' cp_mrf2d(Z_potts, mrfi(1), theta_potts*0.1, c(57,31))
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \url{https://arxiv.org/abs/2006.00383}
#'
#' @export
cp_mrf2d <- function(Z, mrfi, theta, pos){
  C <- max(Z, na.rm = TRUE)
  N <- nrow(Z); M <- ncol(Z)
  if(!any(is.na(Z))){
    probs <- conditional_probabilities_mrf(Z, pos, mrfi@Rmat,
                                           theta, N, M, nrow(mrfi@Rmat), C)
  } else {
    sub_mat <- !is.na(Z)
    probs <- conditional_probabilities_mrf_sub(Z, sub_mat, pos, mrfi@Rmat,
                                           theta, N, M, nrow(mrfi@Rmat), C)
  }
  names(probs) <- 0:C
  return(probs)
}
