#' @name pl_mrf2d
#' @author Victor Freguglia
#' @title Pseudo-likelihood function for MRFs on 2d lattices.
#'
#' @description Computes the pseudo-likelihood function of a Markov Random Field
#'  on a 2-dimensional lattice.
#'
#' @param Z A `matrix` with integers in {0,...,C}.
#' @param mrfi A \code{\link[=mrfi-class]{mrfi}} object representing the
#'  interaction structure.
#' @param theta A 3-dimensional array describing potentials. Slices represent
#'  interacting positions, rows represent pixel values and columns represent
#'  neighbor values. As an example: `theta[1,3,2]` has the potential pairs of
#'  values 0,2 in the second relative position of `mrfi`.
#' @param log_scale A `logical` value indicatin g whether the returned value
#'  should be in logarithmic scale.
#'
#' @return A `numeric` with the pseudo-likelihood value.
#'
#' @details The pseudo-likelihood function is defined as the product of
#' conditional distributions:
#'
#' \deqn{\prod_i P(Z_i | Z_{{N}_i}).}
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


#' @name fit_pl
#' @author Victor Freguglia
#' @title Maximum Pseudo-likelihood fitting of MRFs on 2d lattices.
#'
#' @description Parameter estimation for Markov random fields via
#' Pseudo-Likelihood function optimization. See
#' \code{\link[=pl_mrf2d]{pl_mrf2d}} for more information on the
#' Pseudo-Likelihood function.
#'
#' @param Z A `matrix` object containing the observed MRF.
#' @param mrfi A \code{\link[=mrfi-class]{mrfi}} object representing the
#'  interaction structure.
#' @param family The family of restrictions to potentials. It can be:
#'  * `"onepar"`: Models with a single parameter. It is a single value for any
#'  equal valued pairs and the same value with opposite signal for different
#'  valued pairs.
#' @param init The initial value to be used in the optimization. It can be:
#'  * An `array` with vaid potential values according to `family`.
#'  * A `numeric` of length 1 if `family` is `onepar`.
#'  * `0`. If set to `0` an array with `0`` in all entries is created.
#' @param optim_args Additional parameters passed in `optim()` function call.
#' @param return_optim `logical` indicating whether information from the
#' `optim()` call are returned.
#' @return A `list` object with elements:
#'  * `theta`: The array of estimated potential values.
#'  * `value`: The optimal pseudo-likelihood value.
#'  * `opt.xxx`(if `return_optim` is `TRUE`): Information returned by the
#'   `optim()` function used for the optimization.
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
  C <- length(unique(as.vector(Z))) - 1
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
                             control = list(fnscale = -1)),
                        optim_args))
  out <- list(theta = vec_to_array(o$par, family, C, n_R),
              value = o$value)
  if(return_optim) {out <- c(out, opt = o)}
  return(out)
}
