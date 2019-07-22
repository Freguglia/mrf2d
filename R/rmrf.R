#' @name rmrf2d
#' @author Victor Freguglia
#' @title Sampling of Markov Random Fields on 2d lattices
#'
#' @description Performs pixelwise updates based on conditional distributions
#'  to sample from a Markov random field.
#'
#' @inheritParams pl_mrf2d
#' @param init_Z One of two options:
#'  * A `matrix` object with the initial field configuration. Its
#'  valuesmust be integers in `0,...,C`.
#'  * A length 2 `numeric` vector with the lattice dimensions.
#' @param cycles The number of complete (all pixels) updates to be done.
#' @param method Method used to perform pixel-wise updates.
#'  * `'gibbs'`: Gibbs Sampler method. Samples from conditional distribution.
#'
#' @return A `matrix` with the sampled field.
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
#'
#' @note As in any Gibbs Sampling scheme, a large number of cycles may be
#'  required to achieve the target distribution, specially for strong
#'  interaction systems.
#'
#' @examples
#' # Sample using specified lattice dimension
#' Z <- rmrf2d(c(150,150), mrfi(1), theta_potts)
#'
#' #Sample using itial configuration
#' Z2 <- rmrf2d(Z, mrfi(1), theta_potts)
#'
#' # View results
#' dplot(Z)
#' dplot(Z2)
#'
#' @export
rmrf2d <- function(init_Z, mrfi, theta, cycles = 60, method = "gibbs"){
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
