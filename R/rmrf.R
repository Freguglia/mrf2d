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
#' @param sub_region `NULL` if the whole lattice is considered or a `logical`
#' `matrix` with `TRUE` for pixels in the considered region.
#' @param fixed_region `NULL` if the whole lattice is to be sampled or a
#' `logical` `matrix` with `TRUE` for pixels to be considered fixed. Fixed
#' pixels are not updated in the Gibbs Sampler.
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
#'  A MRF can be sampled in a non-rectangular region of the lattice with the use of
#'  the `sub_region` argument or by setting pixels to `NA` in the initial
#'  configuration `init_Z`. Pixels with `NA` values in `init_Z` are completely
#'  disconsidered from the conditional probabilities and have the same effect as
#'  setting `sub_region = is.na(init_Z)`. If `init_Z` has `NA` values,
#'  `sub_region` is ignored and a warning is produced.
#'
#'  A specific region can be kept constant during the Gibbs Sampler by using the
#'  `fixed_region` argument. Keeping a subset of pixels constant is useful when
#'  you want to sample in a specific region of the image conditional to the
#'  rest, for example, in texture synthesis problems.
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
#' \donttest{
#' Z2 <- rmrf2d(Z, mrfi(1), theta_potts)
#'
#' # View results
#' dplot(Z)
#' dplot(Z2)
#'
#' # Using sub-regions
#' subreg <- matrix(TRUE, 150, 150)
#' subreg <- abs(row(subreg) - 75) + abs(col(subreg) - 75) <= 80
#' # view the sub-region
#' dplot(subreg)
#'
#' Z3 <- rmrf2d(c(150,150), mrfi(1), theta_potts, sub_region = subreg)
#' dplot(Z3)
#'
#' # Using fixed regions
#' fixreg <- matrix(as.logical(diag(150)), 150, 150)
#' # Set initial configuration: diagonal values are 0.
#' init_Z4 <- Z
#' init_Z4[fixreg] <- 0
#'
#' Z4 <- rmrf2d(init_Z4, mrfi(1), theta_potts, fixed_region = fixreg)
#' dplot(Z4)
#'
#' # Combine fixed regions and sub-regions
#' Z5 <- rmrf2d(init_Z4, mrfi(1), theta_potts,
#' fixed_region = fixreg, sub_region = subreg)
#' dplot(Z5)
#' }
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \doi{10.18637/jss.v101.i08}.
#'
#' \code{\link[=rmrf2d_mc]{rmrf2d_mc}} for generating multiple points of a
#' Markov Chain to be used in Monte-Carlo methods.
#'
#' @export
rmrf2d <- function(init_Z, mrfi, theta, cycles = 60, sub_region = NULL, fixed_region = NULL){
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

  if(is.null(sub_region)){ # sub_region is NULL
    if(any(is.na(init_Z))){
      sub_region <- !is.na(init_Z)
    }
    ###########################################################
  } else if(is.matrix(sub_region)) { #sub_region is matrix
    if(!identical(dim(init_Z), dim(sub_region))){
      stop("'init_Z' and 'sub_region' must have the same dimension.")
    }
    if(any(is.na(init_Z))){ # and there are NAs
      warning("'init_Z' has NA values and 'sub_region' was defined. Using non-NA values in 'init_Z' as sub-region and ignoring 'sub_region'")
      sub_region <- !is.na(init_Z)
    }
    #########################################################
  } else { # sub_region is wrong
    stop("'sub_region' must be either NULL or a logical matrix.")
  }

  if(is.matrix(fixed_region)) { #fixed_region is matrix
    if(!identical(dim(init_Z), dim(fixed_region))){
      stop("'init_Z' and 'fixed_region' must have the same dimension.")
    }
    #########################################################
  } else if(!is.null(fixed_region)) { # sub_region is wrong
    stop("'fixed_region' must be either NULL or a logical matrix.")
  }

  if(is.null(fixed_region) && !is.null(sub_region)){
    fixed_region <- matrix(FALSE, nrow(init_Z), ncol(init_Z))
  }

  if(!is.null(fixed_region) && is.null(sub_region)){
    sub_region <- matrix(TRUE, nrow(init_Z), ncol(init_Z))
  }

  if(!is.null(fixed_region)){
    if(any(fixed_region[!sub_region])){
      warning("Some pixels in the 'fixed_region' are not part of the 'sub_region', they will be ignored.")
    }
  }

  theta <- sanitize_theta(theta)

  R <- mrfi@Rmat

  null_interactions <- apply(theta, MARGIN = 3, function(m) all(m == 0))
  theta <- theta[,,!null_interactions]
  R <- R[!null_interactions,]

  if(is.null(sub_region) && is.null(fixed_region)){
    return(gibbs_sampler_mrf2d(init_Z, R, theta, cycles))
  } else{
    ret <- gibbs_sampler_mrf2d_sub(init_Z, sub_region, fixed_region, R, theta, cycles)
    ret[!sub_region] <- NA
    return(ret)
  }

}

#' @name rmrf2d_mc
#' @author Victor Freguglia
#' @title Markov Chain sampling of MRFs for Monte-Carlo methods
#'
#' @description Generates a Markov Chain of random fields and returns the
#' sufficient statistics for each of the observations.
#'
#' This function automatizes the process of generating a random sample of MRFs
#' to be used in Monte-Carlo methods by wrapping \code{\link[=rmrf2d]{rmrf2d}}
#' and executing it multiple time while storing sufficient statistics instead
#' of the entire lattice.
#'
#' @inheritParams rmrf2d
#' @inheritParams smr_stat
#' @param burnin Number of cycles iterated before start collecting sufficient
#' statistics.
#' @param cycles Number of cycles between collected samples.
#' @param nmc Number of samples to be stored.
#' @param verbose `logical` indicating whether to print iteration number or not.
#'
#' @note Fixed regions and incomplete lattices are not supported.
#'
#' @return A matrix where each row contains the vector of sufficient statistics
#' for an observation.
#'
#' @examples
#' rmrf2d_mc(c(80, 80), mrfi(1), theta_potts, family = "oneeach", nmc = 8)
#'
#' @export
rmrf2d_mc <- function(init_Z, mrfi, theta, family,
                      nmc = 100, burnin = 100, cycles = 4,
                      verbose = interactive()){

  # Initialize samples and matrix to store results
  zt <- rmrf2d(init_Z, mrfi, theta, cycles = burnin)
  smrs <- smr_stat(zt, mrfi, family)
  M <- matrix(numeric(), nrow = nmc, ncol = length(smrs))

  # Generate Markov Chain of Random Fields
  for(i in seq_len(nmc)){
    zt <- rmrf2d(zt, mrfi, theta, cycles = cycles)
    # Store sufficient stastics
    M[i,] <- smr_stat(zt, mrfi, family)
    if(verbose) cat("\r", i)
  }
  if(verbose) cat("\r", "Done!", "\n")

  return(M)
}
