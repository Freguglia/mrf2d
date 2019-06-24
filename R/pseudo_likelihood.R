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

fitpl_onepar_mrf2d <- function(Z, mrfi, init = 0, return_optim = FALSE,
                               optim_args = list(method = "BFGS")){
  if(!length(init) == 1 | !is.numeric(init)) {
    stop("Argument 'init' must be a length 1 numeric.")
  }
  C <- length(unique(as.vector(Z))) - 1
  R <- mrfi@Rmat
  n_R <- nrow(R)
  pl_value <- function(par){
    # Create theta array based on one parameter
    theta <- array(diag(C+1)*par*2 - par, dim = c(C+1,C+1,n_R))
    # Compute the pseudo-likelihood value
    return(pl_mrf2d(Z, mrfi, theta))
  }
  o <- do.call(optim, c(list(par = init,
                             fn = pl_value,
                             control = list(fnscale = -1)),
                        optim_args))
  out <- list(theta = array(diag(C+1)*o$par*2 - o$par, dim = c(C+1,C+1,n_R)),
              value = o$value)
  if(return_optim) {out <- c(out, o)}
  return(out)
}
