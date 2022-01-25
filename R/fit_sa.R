#' @name fit_sa
#' @author Victor Freguglia
#' @title Stochastic Approximation fitting of MRFs on 2d lattices
#'
#' @description Estimates the parameters of a MRF by successively sampling from
#'  a parameter configuration and updating it by comparing the sufficient statistics
#'  of the sampled field and the observed field.
#'
#'  This method aims to find the parameter value where the gradient of the
#'  likelihood function is equal to zero.
#'
#'
#' @inheritParams fit_pl
#' @inheritParams rmrf2d
#' @param gamma_seq A `numeric` vector with the sequence of constants
#' used in each step \eqn{\gamma_t}.
#' @param refresh_each An integer with the number of iterations taken before a
#' complete refresh (restart from a random state). This prevents the sample from
#' being stuck in a mode for too long. Defaults to `length(gamma_seq) + 1` (no
#' refresh happens).
#' @param refresh_cycles An integer indicating how many Gibbs Sampler cycles are
#' performed when a refresh happens. Larger is usually better, but slower.
#' @param verbose `logical` indicating whether the iteration number is printed
#' during execution.
#'
#' @return A `mrfout` object with the following elements:
#'  * `theta`: The estimated `array` of potentials.
#'  * `mrfi`: The interaction structure considered.
#'  * `family`: The parameter restriction family considered.
#'  * `method`: The estimation method (`"Stochastic Approximation"`).
#'  * `metrics`: A `data.frame` containing the the euclidean distance between
#'  the sufficient statics computed for `Z` and the current sample.
#'
#' @details
#'
#' The stochastic approximation method consists of, given an observed field `Z`,
#' and a starting parameters configuration \eqn{\theta_0}, successively sample
#' a field \eqn{Z_t} from the current parameter configuration and estimate the
#' direction of the  gradient of the likelihood function by comparing the
#' sufficient statistics in the current sample and the observed field.
#'
#' The solution is updated by moving in the estimated direction with a predefined
#' step size \eqn{\gamma_t}, a new field \eqn{Z_{t+1}} is sampled using the new
#' parameter configuration and \eqn{Z_t} as an initial value, and the process is
#' repeated.
#'
#' \deqn{\theta_{t+1} = \theta_t - \gamma_t(T(Z_t) - T(Z)),}
#' where \eqn{T(Z)} is the sufficient statistics for the reference field,
#' \eqn{T(Z_t)} is the sufficient statistics for a field sampled from
#' \eqn{\theta_t}.
#'
#' `gamma_seq` is normalized internally by diving values by `length(Z)`, so the
#' choice of the sequence is invariant to the lattice dimensions. Typically, a
#' sequence like `seq(from = 1, to = 0, length.out = 1000)` should be used for
#' defining a sequence with `1000` steps. Some tuning of this sequence is
#' required.
#'
#' @note
#'   Stochastic Approximation is called "Controllable Simulated Annealing" in
#' some references.
#'
#' Examples where Stochastic Approximation is used with MRFs are
#' \insertCite{gimel_sa}{mrf2d}, \insertCite{CR}{mrf2d}.
#'
#' @references
#' \insertRef{wiki_sa}{mrf2d}
#'
#' \insertAllCited{}
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' fit1 <- fit_sa(Z_potts, mrfi(1), family = "oneeach", gamma_seq = seq(1, 0, length.out = 50))
#' # Estimated parameters
#' fit1$theta
#' # A visualization of estimated gradient norm over iterations.
#' plot(fit1$metrics, type = "l")
#'
#' fit_sa(Z_potts, mrfi(1), family = "oneeach", gamma_seq = seq(1, 0, length.out = 50))
#' }
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \doi{10.18637/jss.v101.i08}.
#'
#' @export
fit_sa <- function(Z, mrfi, family = "onepar", gamma_seq, init = 0, cycles = 5,
                   refresh_each = length(gamma_seq)+1, refresh_cycles = 60,
                   verbose = interactive()){

  if(!family %in% mrf2d_families){
    stop("'", family, "' is not an implemented family.")
  }
  if(!is.numeric(init)) {
    stop("Argument 'init' must be numeric.")
  }

  C <- length(na.omit(unique(as.vector(Z)))) - 1
  R <- mrfi@Rmat
  n_R <- nrow(R)
  gamma_seq <- gamma_seq/prod(dim(Z))

  if(is.vector(init)) {
    if(identical(init, 0)){
      init <- array(0, dim = c(C+1, C+1, n_R))
    } else {
      init <- vec_to_array(init, family, C, n_R)
    }
  } else if(!is_valid_array(init, family)) {
    stop("'init' array is incompatible with family '", family,"'")
  }

  arr_Z <- table_relative_3d(Z, mrfi@Rmat, C)
  S <- suf_stat(arr_Z, family)
  d <- numeric(length(gamma_seq))
  Zseq <- thetaseq <- matrix(NA, nrow = length(gamma_seq), ncol = length(S))
  is_sub <- any(is.na(Z))
  if(is_sub){
    subr <- !is.na(Z)
  }
  # Initialize
  theta_t <- array_to_vec(init, family)
  if(!is_sub){
    Z_t <- rmrf2d(dim(Z), mrfi, vec_to_array(theta_t, family, C, n_R), cycles)
  } else {
    Z_t <-  rmrf2d(dim(Z), mrfi, vec_to_array(theta_t, family, C, n_R), cycles,
                   sub_region = subr)
  }
  # Iterate
  for(t in seq_along(gamma_seq)){
    cat(ifelse(verbose, paste("\r Iteration:", t), ""))
    arr_Z_t <- table_relative_3d(Z_t, mrfi@Rmat, C)
    S_t <- suf_stat(arr_Z_t, family)
    theta_t <- theta_t - gamma_seq[t]*(S_t - S)
    if(t%%refresh_each == 0){
      if(!is_sub){
        Z_t <- rmrf2d(dim(Z), mrfi, vec_to_array(theta_t, family, C, n_R),
                      refresh_cycles)
      } else {
        Z_t <-  rmrf2d(dim(Z), mrfi, vec_to_array(theta_t, family, C, n_R),
                       refresh_cycles, sub_region = subr)
      }
    } else {
      Z_t <- rmrf2d(Z_t, mrfi, vec_to_array(theta_t, family, C, n_R), cycles)
    }
    Zseq[t,] <- S_t
    thetaseq[t,] <- theta_t
    d[t] <- sqrt(sum((S_t - S)^2))
  }
  cat(ifelse(verbose, "\n", ""))
  theta_out <- vec_to_array(theta_t, family, C, n_R)
  dimnames(theta_out)[[3]] <- mrfi_to_char(mrfi)
  out <- list(theta = theta_out,
              mrfi = mrfi,
              family = family,
              method = "Stochastic Approximation",
              metrics = data.frame(t = seq_along(gamma_seq), distance = d),
              Zseq = Zseq,
              thetaseq = thetaseq,
              Z = Z,
              ncycles = length(gamma_seq))
  class(out) <- "mrfout"
  return(out)
}
