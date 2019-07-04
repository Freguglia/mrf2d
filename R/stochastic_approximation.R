#' @name fit_sa
#' @author Victor Freguglia
#' @title Stochastic Approximation algorithm for MRFs on 2d lattices
#'
#' @description Parameter estimation for Markov Random Fieds via Stochastic
#'  Approximation algorithm.
#'
#' @inheritParams fit_pl
#' @param gamma_seq A `numeric` vector with the sequence of constants \eqn{\gamma_t}
#' @param steps Number of Gibbs Sampler updates between samples.
#'
#' @return A `list` object with the following elements:
#'  * `theta`: The estimated `array` of potentials.
#'  * `metrics`: A `data.frame` containing the the euclidean distance between
#'  the sufficient statics computed for `Z` and the current sample.
#'
#' @export
fit_sa <- function(Z, mrfi, family = "onepar", gamma_seq, init = 0, steps = 5){

  if(!family %in% mrf2d_families){
    stop("'", family, "' is not an implemented family.")
  }
  if(!is.numeric(init)) {
    stop("Argument 'init' must be numeric.")
  }

  C <- length(unique(as.vector(Z))) - 1
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

  arr_Z <- table_relative_3d(Z, mrfi@Rmat, C, FALSE)
  S <- suf_stat(arr_Z, family)
  d <- numeric(length(gamma_seq))
  # Initialize
  theta_t <- array_to_vec(init, family)
  Z_t <- rmrf2d(dim(Z), mrfi, vec_to_array(theta_t, family, C, n_R), steps)
  # Iterate
  for(t in seq_along(gamma_seq)){
    arr_Z_t <- table_relative_3d(Z_t, mrfi@Rmat, C, FALSE)
    S_t <- suf_stat(arr_Z_t, family)
    theta_t <- theta_t - gamma_seq[t]*(S_t - S)
    Z_t <- rmrf2d(Z_t, mrfi, vec_to_array(theta_t, family, C, n_R), steps)
    d[t] <- sqrt(sum((S_t - S)^2))
  }
  return(list(theta = vec_to_array(theta_t, family, C, n_R),
              metrics = data.frame(t = seq_along(gamma_seq), distance = d)))
}
