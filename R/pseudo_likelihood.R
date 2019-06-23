#' @name pl_mrf2d
#' @author Victor Freguglia
#' @title Pseudo-likelihood function for MRFs on 2d lattices.
#'
#' @description Computes the pseudo-likelihood function of a Markov Random Field
#'  on a 2-dimensional lattice.
#'
#' @details Add description of Pseudo-Likelihood. Product of conditionals.
#'
#' @param Z A `matrix` with integers in {0,...,C}.
#' @param mrfi A \code{\link[=mrfi-class]{mrfi}} object representing the
#'  interaction structure.
#' @param theta A 3-dimensional array describing potentials. Slices represent
#'  interacting positions, rows represent pixel values and columns represent
#'  neighbor values. As an example: `theta[1,3,2]` has the potential pairs of
#'  values 0,2 in the second relative position of `mrfi`.
#' @param log_scale A `logical` value indicating whether the returned value
#'  should be in logarithmic scale.
#'
#' @return A `numeric` with the pseudo-likelihood value.
#'
#' @export
pl_mrf2d <- function(Z, mrfi, theta, log_scale = TRUE){
  R <- mrfi@Rmat
  pl <- log_pl_mrf(Z, R, theta)
  if(log_scale) {return(pl)
  } else {
    return(exp(pl))
  }
}
