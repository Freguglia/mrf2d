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
#'  valuesmust be integers in `{0,...,C}`.
#'  * A length 2 `numeric` vector with the lattice dimensions.
#' @param cycles The number of updates to be done (for each each pixel).
#' @param sub_lattice `NULL` if the whole lattice is considered or a `logical`
#' `matrix` with `TRUE` for pixels in the considered region.
#' @param mask_na `logical` indicating whether pixels not in the sub-lattice
#' should be set to NA.
#'
#' @return A `matrix` with the sampled field.
#'
#' @details This function implements a Gibbs Sampling scheme to sample from
#' a Markov random field by iteratively sampling pixel values from the
#' conditional distribution
#'  \deqn{P(Z_i | Z_{{N}_i}, \theta).}
#'
#'  A cycle means exactly one update to each pixel. The order pixels are
#'  sampled is randomized within each cycle.
#'
#'  If `init_Z` is passed as a length 2 vector with lattice dimensions, the
#'  initial field is sampled from independent discrete uniform distributions in
#'  `{0,...,C}`. The value of C is obtained from the number of rows/columns of
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
#' # Using sub-lattices
#' sublat <- matrix(TRUE, 150, 150)
#' sublat <- abs(row(sublat) - 75) + abs(col(sublat) - 75) <= 80
#' # view the sub-lattice region
#' dplot(sublat)
#'
#' Z3 <- rmrf2d(c(150,150), mrfi(1), theta_potts, sub_lattice = sublat)
#' dplot(Z3)
#'
#' @export
rmrf2d <- function(init_Z, mrfi, theta, cycles = 60, sub_lattice = NULL, mask_na = TRUE){
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

  if(is.null(sub_lattice)){ # sub_lattice is NULL
    if(any(is.na(init_Z))){
      sub_lattice <- !is.na(init_Z)
    }
    ###########################################################
  } else if(is.matrix(sub_lattice)) { #sub_lattice is matrix
    if(!identical(dim(init_Z), dim(sub_lattice))){
      stop("'init_Z' and 'sub_lattice' must have the same dimension.")
    }
    if(any(is.na(init_Z))){ # and there are NAs
      warning("'init_Z' has NA values and 'sub_lattice' was defined. Using non-NA values in 'init_Z' as sub-lattice and ignoring 'sub_lattice'")
      sub_lattice <- !is.na(init_Z)
    }
    #########################################################
  } else { # sub_lattice is wrong
    stop("'sub_lattice' must be either NULL or a logical matrix.")
  }

  theta <- sanitize_theta(theta)

  R <- mrfi@Rmat

  null_interactions <- apply(theta, MARGIN = 3, function(m) all(m == 0))
  theta <- theta[,,!null_interactions]
  R <- R[!null_interactions,]

  if(is.null(sub_lattice)){
    return(gibbs_sampler_mrf2d(init_Z, R, theta, cycles))
  } else{
    ret <- gibbs_sampler_mrf2d_sub(init_Z, sub_lattice, R, theta, cycles)
    if(mask_na){
      ret <- ifelse(!sub_lattice, NA, ret)
    }
    return(ret)
  }


}
