#' @name rmrf2d
#' @title Sampling of Markov Random Fields on 2d lattices
#'
#' @description Performs pixel-wise updates based on conditional distributions
#'  to sample from a Markov random field. The order pixels are updated is
#'  randomized in each step.
#'
#' @param init_Z A `matrix` object with the initial field configuration. Its
#'  valuesmust be integers in 0,...,C. A length 2 `numeric` vector with lattice
#'  dimensions can be used to start from a random configuration.
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
rmrf2d <- function(init_Z, mrfi, theta, steps = 10, method = "gibbs"){
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
    return(gibbs_sampler_mrf2d(init_Z, R, theta, steps))
  } else {
    valid_methods <- c("'gibbs'")
    stop("Invalid 'method' argument. Current support methods are: ",
         paste(valid_methods, collapse = " "))
  }

}

#' @export
prmrf2d <- function(init_Z, mrfi, theta, steps = 10, method = "gibbs", nblocks = 2){
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
  Z <- init_Z

  len <- length(Z)

  blocks <- sort(rep_len(1:(2*nblocks), len))
  bl_lis <- split((1:len) -1, blocks)
  blmat <- matrix(blocks, nrow(Z), ncol(Z), byrow = T)
  odd_bl <- seq(1, 2*nblocks-1, 2)
  even_bl <- seq(2, 2*nblocks, 2)

  if(method == "gibbs"){
    for(i in seq_len(steps)){
      #Process odd blocks
      odd_res <- foreach(bl = odd_bl) %dopar% {
        pgibbs_sampler_mrf2d(Z, R, theta, 1, subset = bl_lis[[bl]])
      }

      # Update Z with the new odd blocks
      for(j in seq_along(odd_bl)){
        idx <- blmat == odd_bl[j]
        Z[idx] <- odd_res[[j]][idx]
      }

      #Process even blocks
      odd_res <- foreach(bl = even_bl) %dopar% {
        pgibbs_sampler_mrf2d(Z, R, theta, 1, subset = bl_lis[[bl]])
      }

      # Update Z with the new odd blocks
      for(j in seq_along(even_bl)){
        idx <- blmat == even_bl[j]
        Z[idx] <- odd_res[[j]][idx]
      }

    }
    return(Z)

  } else {
    valid_methods <- c("'gibbs'")
    stop("Invalid 'method' argument. Current support methods are: ",
         paste(valid_methods, collapse = " "))
  }

}

# benchmark(sample_in = rmrf2d(c(100,100), mrfi(), theta_potts, steps = 6),
#           sample_out = prmrf2d(c(100,100), mrfi(), theta_potts, steps = 6),
#           replications = 10)
