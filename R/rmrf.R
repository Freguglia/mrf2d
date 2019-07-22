#' @name rmrf2d
#' @author Victor Freguglia
#' @title Sampling of Markov Random Fields on 2d lattices
#'
#' @description Performs pixelwise updates based on conditional distributions
#'  to sample from a Markov random field.
#'
#'
#' @param init_Z One of two options:
#'  * A `matrix` object with the initial field configuration. Its
#'  valuesmust be integers in `0,...,C`.
#'  * A length 2 `numeric` vector with the lattice dimensions.
#' @param mrfi A \code{\link[=mrfi-class]{mrfi}} object representing the
#'  interaction structure.
#' @param theta A 3-dimensional array describing potentials. Slices represent
#'  interacting positions, rows represent pixel values and columns represent
#'  neighbor values. As an example: `theta[1,3,2]` has the potential pairs of
#'  values 0,2 in the second relative position of `mrfi`.
#' @param cycles The number of complete (all pixels) updates to be done.
#' @param method Method used to perform pixel-wise updates.
#'  * `'gibbs'`: Gibbs Sampler method. Samples from conditional distribution.
#'
#' @return
#'
#' @details This function implements a Gibbs Sampling scheme to sample from
#' a Markov random field by iteratively sampling pixel values from the
#' conditional distribution
#'  \deqn{P(Z_i | Z_{{N}_i}).}
#'
#'  A cycle means exactly one update to each pixel. The order pixels are
#'  sampled is randomized within each cycle.
#'
#'  If `init_Z` is passed as a length 2 vector with lattice dimensions, the
#'  initial field is sampled from independent discrete uniform distributions in
#'  `0,...,C`. The value of C is obtained from the number of rows/columns of
#'  `theta`.
#' @export
rmrf2d <- function(init_Z, mrfi, theta, cycles = 10, method = "gibbs"){
  # Check validity of the input
  if(!is.matrix(init_Z)) {
    if(is.numeric(init_Z) & is.vector(init_Z)) {
      if(length(init_Z) == 2 & min(init_Z) > 0) {
        .space <- 0:(dim(theta)[1] - 1)
        init_Z <- matrix(sample(.space, prod(init_Z),replace = TRUE),
                         nrow = init_Z[1], ncol = init_Z[2])
      } else {
        stop("Argument 'init_Z' must be either a matrix or a length 2 numeric vector with lattice dimensions.")
      }
    } else {
      stop("Argument 'init_Z' must be either a matrix or a length 2 numeric vector with lattice dimensions.")
    }
  }
  theta <- sanitize_theta(theta)

  R <- mrfi@Rmat

  if(method == "gibbs"){
    return(gibbs_sampler_mrf2d(init_Z, R, theta, cycles))
  } else {
    valid_methods <- c("'gibbs'")
    stop("Invalid 'method' argument. Current support methods are: ",
         paste(valid_methods, collapse = " "))
  }

}
