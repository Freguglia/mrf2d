#' @name pl_mrf2d
#' @author Victor Freguglia
#' @title Pseudo-likelihood function for MRFs on 2d lattices
#'
#' @description Computes the pseudo-likelihood function of a Markov Random Field
#'  on a 2-dimensional lattice.
#'
#' @param Z A `matrix` with integers in `{0,...,C}`.
#' @param mrfi A \code{\link[=mrfi-class]{mrfi}} object representing the
#'  interaction structure.
#' @param theta A 3-dimensional array describing potentials. Slices represent
#'  interacting positions, rows represent pixel values and columns represent
#'  neighbor values. As an example: `theta[1,3,2]` has the potential for the
#'  pair of values 0,2 observed in the second relative position of `mrfi`.
#' @param log_scale A `logical` value indicating whether the returned value
#'  should be in logarithmic scale.
#'
#' @return A `numeric` with the pseudo-likelihood value.
#'
#' @details The pseudo-likelihood function is defined as the product of
#' conditional probabilities of each pixel given its neighbors:
#'
#' \deqn{\prod_i P(Z_i | Z_{{N}_i}, \theta).}
#'
#' @examples
#' pl_mrf2d(Z_potts, mrfi(1), theta_potts)
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \doi{10.18637/jss.v101.i08}.
#'
#' @export
pl_mrf2d <- function(Z, mrfi, theta, log_scale = TRUE){
  if(any(is.na(Z))){
    return(pl_sub(Z, mrfi, theta, log_scale))
  } else {
    return(pl_nosub(Z, mrfi, theta, log_scale))
  }
}

pl_nosub <- function(Z, mrfi, theta, log_scale){
  R <- mrfi@Rmat
  pl <- log_pl_mrf(Z, R, theta)
  if(log_scale) {return(pl)
  } else {
    return(exp(pl))
  }
}

pl_sub <- function(Z, mrfi, theta, log_scale){
  sub_mat <- !is.na(Z)
  R <- mrfi@Rmat
  pl <- log_pl_mrf_sub(Z, sub_mat, R, theta)
  if(log_scale) {return(pl)
  } else {
    return(exp(pl))
  }
}


#' @name fit_pl
#' @author Victor Freguglia
#' @title Maximum Pseudo-likelihood fitting of MRFs on 2d lattices.
#'
#' @description Parameter estimation for Markov random fields via
#' Pseudo-Likelihood function optimization. See
#' \code{\link[=pl_mrf2d]{pl_mrf2d}} for information on the
#' Pseudo-Likelihood function.
#'
#' @inheritParams pl_mrf2d
#' @param Z A `matrix` object containing the observed MRF. `NA` values can be
#' used to create a subregion of the lattice for non-rectangular data.
#' @param family The family of parameter restrictions to potentials. Families
#' are:
#'   `'onepar'`, `'oneeach'`, `'absdif'`, `'dif'` or `'free'`.
#' See \code{\link[=mrf2d-family]{mrf2d-familiy}}.
#' @param init The initial value to be used in the optimization. It can be:
#'  * A valid `array` of parameter values according to `family`.
#'  * `0`. If set to `0` an array with `0`` in all entries is created.
#' @param optim_args Additional parameters passed to `optim()`.
#' @param return_optim `logical` indicating whether information from the
#' `optim()` call are returned.
#'
#' @return An object of class `mrfout` with elements:
#'  * `theta`: The estimated array of potential values.
#'  * `mrfi`: The interaction structure considered.
#'  * `family`: The parameter restriction family considered.
#'  * `method`: The estimation method (`"Pseudolikelihood"`).
#'  * `value`: The optimal pseudo-likelihood value.
#'  * `opt.xxx`(if `return_optim` is `TRUE`): Information returned by the
#'   `optim()` function used for the optimization.
#'
#'
#' @examples
#' fit_pl(Z_potts, mrfi(1), family = "onepar")
#' fit_pl(Z_potts, mrfi(1), family = "oneeach")
#' fit_pl(Z_potts, mrfi(2), family = "onepar")
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \doi{10.18637/jss.v101.i08}.
#'
#' @importFrom stats optim
#' @export
fit_pl <- function(Z, mrfi, family = "onepar", init = 0,
                   optim_args = list(method = "BFGS"),
                   return_optim = FALSE){

  if(!family %in% mrf2d_families){
    stop("'", family, "' is not an implemented family.")
  }
  if(!is.numeric(init)) {
    stop("Argument 'init' must be numeric.")
  }
  C <- length(na.omit(unique(as.vector(Z)))) - 1
  R <- mrfi@Rmat
  n_R <- nrow(R)

  if(is.vector(init)) {
    if(identical(init, 0)){
      init <- array(0, dim = c(C+1, C+1, n_R))
    } else {
      init <- vec_to_array(init, family, C, n_R)
    }
  } else if(!is_valid_array(init, family)) {
    stop("'init' array is incompatible with family '", family,"'")
  }
  init <- array_to_vec(init, family)
  pl_value <- function(par){
    # Create theta array based on one parameter
    theta <- vec_to_array(par, family, C, n_R)
    # Compute the pseudo-likelihood value
    return(pl_mrf2d(Z, mrfi, theta))
  }
  o <- do.call(optim, c(list(par = init,
                             fn = pl_value,
                             control = list(fnscale = -1),
                             hessian = return_optim),
                        optim_args))
  theta_out = vec_to_array(o$par, family, C, n_R)
  dimnames(theta_out)[[3]] <- mrfi_to_char(mrfi)
  out <- list(theta = theta_out,
              mrfi = mrfi,
              family = family,
              method = "Pseudolikelihood",
              value = o$value,
              Z = Z)
  if(return_optim) {out <- c(out, opt = o)}
  class(out) <- "mrfout"
  return(out)
}

# Gradient for Pseudolikelihood under 'free' family
grad_pl_free_nosub <- function(theta_vec, z, mrfi){
  ns <- cohist(z, mrfi)
  C <- dim(ns)[2] - 1
  theta_arr <- expand_array(theta_vec, "free", mrfi, C)
  wprob <- gradient_crossed_free(z, mrfi@Rmat, theta_arr)
  out <- 2*ns - wprob
  out[1,1,] <- 0
  return(out)
}
