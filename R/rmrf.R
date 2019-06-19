#' @name rmrf
#' @title Sampling of Markov Random Fields
#'
#' @description Performs pixel-wise updates based on conditional distributions
#'  to sample from a Markov random field. The order pixels are updated is
#'  randomized in each step.
#'
#' @param init_Z A `matrix` object with the initial field configuration. Its values
#'  must be integers in 0,...,C.
#' @param mrfi A \code{\link[=mrfi-class]{mrfi}} object representing the
#'  interaction structure.
#' @param theta A 3-dimensional array describing potentials. Slices represent
#'  interacting positions, rows represent pixel values and columns represent
#'  neighbor values. As an example: `theta[1,3,2]` has the potential pairs of
#'  values 0,2 in the second relative position of `mrfi`.
#' @param steps The number of complete (all pixels) updates to be done.
#' @param method Method used to perform pixel-wise updates.
#'  * `'gibbs'`: Gibbs Sampler method. Samples from conditional distribution.
#' @export
rmrf <- function(init_Z, mrfi, theta, steps = 10, method = "gibbs"){
  R <- mrfi@Rmat

  if(method == "gibbs"){
    return(gibbs_sampler_mrf(init_Z, R, theta, steps))
  }
  else {
    valid_methods <- c("'gibbs'")
    stop("Invalid 'method' argument. Current support methods are: ",
         paste(valid_methods, collapse = " "))
  }

}
