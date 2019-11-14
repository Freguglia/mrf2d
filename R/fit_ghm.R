#' @name fit_ghm
#' @author Victor Freguglia
#' @title EM estimation for Gaussian Hidden Markov field
#'
#' @description `fit_ghm` fits a Gaussian Mixture model with hidden components
#' driven by a Markov random field with known parameters. The inclusion of a
#' linear combination of basis functions as a fixed effect is also possible.
#'
#' The algorithm is an implementation of
#' \insertCite{zhang2001segmentation}{mrf2d}.
#'
#'
#' @param Y A matrix of observed (continuous) pixel values.
#' @inheritParams pl_mrf2d
#' @param fixed_fn A list of functions `fn(x,y)` to be considered as a fixed
#'  effect. See \code{\link[=basis_functions]{basis_functions}}.
#' @param equal_vars `logical` indicating if the mixture model has the same
#'  variances in all mixture components.
#' @param init_mus Optional. A `numeric` with length (C+1) with the initial mean
#' estimate for each component.
#' @param init_sigmas Otional. A `numeric` with length (C+1) with initial sample
#' deviation estimate for each component.
#' @param maxiter The maximum number of iterations allowed. Defaults to 100.
#' @param max_dist Defines a stopping condition. The algorithm stops if the
#'  maximum absolute difference between parameters of two consecutive iterations
#'  is less than `max_dist`.
#' @param icm_cycles Number of steps used in the Iterated Conditional Modes
#'  algorithm executed in each interaction. Defaults to 6.
#' @param verbose `logical` indicating wheter to print the progress or not.
#' @param qr The QR decomposition of the design matrix. Used internally.
#'
#' @return A `list` containing:
#'  * `par`: A `data.frame` with \eqn{\mu} and \eqn{\sigma} estimates for each
#' component.
#'  * `fixed`: A `matrix` with the estimated fixed effect in each pixel.
#'  * `Z_pred`: A `matrix` with the predicted component (highest probability) in
#'  each pixel.
#'  * `predicted`: A `matrix` with the fixed effect + the \eqn{\mu} value for
#'  the predicted component in each pixel.
#'  * `iterations`: Number of EM iterations done.
#'
#' @details If either `init_mus` or `init_sigmas` is `NULL` an EM algorithm
#'  considering an independent uniform distriburion for the hidden component is
#'  fitted first and its estimated means and sample deviations are used as
#'  initial values. This is necessary because the algorithm may not converge if
#'  the initial parameter configuration is too far from the maximum likelihood
#'  estimators.
#'
#'  `max_dist` defines a stopping condition. The algorithm will stop if the
#'  maximum absolute difference between (\eqn{\mu} and \eqn{\sigma}) parameters
#'  in consecutive iterations is less than `max_dist`.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Sample a Gaussian mixture with components given by Z_potts
#' # mean values are 0, 1 and 2 and a linear effect on the x-axis.
#' \dontrun{
#' set.seed(2)
#' Y <- Z_potts + rnorm(length(Z_potts), sd = 0.4) +
#'       (row(Z_potts) - mean(row(Z_potts)))*0.01
#' # Check what the data looks like
#' cplot(Y)
#'
#' fixed <- polynomial_2d(c(1,0), dim(Y))
#' fit <- fit_ghm(Y, mrfi = mrfi(1), theta = theta_potts, fixed_fn = fixed)
#' fit$par
#' cplot(fit$fixed)
#' }
#'
#' @importFrom stats lm predict
#' @export
fit_ghm <- function(Y, mrfi, theta, fixed_fn = list(),
                    equal_vars = TRUE,
                    init_mus = NULL,
                    init_sigmas = NULL,
                    maxiter = 100, max_dist = 10^-3,
                    icm_cycles = 6, verbose = TRUE, qr = NULL){

  Rmat <- mrfi@Rmat

  # check if theta's dimensions matches R's
  if(dim(theta)[3] != nrow(Rmat)) {stop("Invalid dimensions.")}
  C <- dim(theta)[1] - 1

  # compute initial fixed effect
  N <- nrow(Y); M <- ncol(Y)
  if(length(fixed_fn) > 0){
    if(is.null(qr)){
      X <- basis_function_df(fixed_fn, N, M,
                             standardize = TRUE)[as.logical(!is.na(Y)),-(1:2)]
      X <- as.matrix(X)
      q <- qr(X)
    } else {
      q <- qr
    }
    S <- Y
    S[!is.na(Y)] <-   qr.fitted(q, Y[!is.na(Y)])
    e <- Y - S
  } else {
    e <- Y
    S <- Y*0
  }

  is_sub <- any(is.na(Y))
  if(is_sub){
    subr <- !is.na(Y)
  }

  # Initialize variables
  if(is.null(init_mus) | is.null(init_sigmas)) {
    mus_old <- seq(min(e, na.rm = TRUE), max(e, na.rm = TRUE), length.out = C+1)
    if(verbose)
      cat("\r Fitting independent mixture to obtain initial parameters. \n")
    ind_fit <- fit_ghm(e, mrfi, theta*0, fixed_fn = fixed_fn, equal_vars,
                       init_mus = seq(min(e, na.rm = TRUE), max(e, na.rm = TRUE),
                                      length.out = C+1),
                       init_sigmas = rep(diff(range(e, na.rm = TRUE))/(2*C), C+1),
                       maxiter, max_dist, icm_cycles, verbose = FALSE)
    mus_old <- ind_fit$par$mu
    sigmas_old <- ind_fit$par$sigma
  } else {
    mus_old <- init_mus
    sigmas_old <- init_sigmas
  }
  Z <- matrix(sample(0:C, size = N*M, replace = TRUE), nrow = nrow(Y))
  Z <- ifelse(is.na(Y), NA, Z)

  # Iterate
  iter <- 0; dist <- Inf
  while(iter < maxiter && dist > max_dist){

    # Compute MAP estimates of Z
    if(!is_sub){
      Z <- icm_gaussian_cpp(e, Rmat, Z, theta, mus_old, sigmas_old, icm_cycles)
      cond_probs <- cprob_ghm_all(Z, Rmat, theta, mus_old, sigmas_old, e)
    } else {
      Z <- icm_gaussian_cpp_sub(e, subr, Rmat, Z, theta, mus_old, sigmas_old, icm_cycles)
      cond_probs <- cprob_ghm_all_sub(Z, subr, Rmat, theta, mus_old, sigmas_old, e)
    }

    # Compute conditional probabilities

    # Update mus, sigmas and S
    ## mus and sigmas
    if(equal_vars){
      mus_new <- apply(cond_probs, MARGIN = 3, function(x)
        sum(x*e, na.rm = TRUE)/sum(x, na.rm = TRUE))
      sigmas_new <- rep(
        sqrt(
          sum(
            simplify2array(
              lapply(mus_new,
                     function(mu) (e - mu)^2))*cond_probs, na.rm = TRUE)/(N*M)),
        C+1)
    }
    else {
      mus_new <- apply(cond_probs, MARGIN = 3, function(x)
        sum(x*e, na.rm = TRUE)/sum(x, na.rm = TRUE))
      sigmas_new <- sapply(1:(C+1), function(l) {
        sum(cond_probs[,,l]*(e - mus_new[l])^2, na.rm = TRUE)/sum(cond_probs[,,l], na.rm = TRUE)
      })
      sigmas_new <- sqrt(sigmas_new)
    }

    ## update S
    if(length(fixed_fn) > 0){
      if(equal_vars){
        mean_est <- apply(cond_probs, MARGIN = c(1,2),
                          function(p_vec){
                            sum(p_vec*mus_new, na.rm = TRUE)
                          })
        pred <- qr.fitted(q, as.vector(Y - mean_est)[!is.na(Y)])
        S <- Y
        S[!is.na(Y)] <- pred
      } else {
        norm_cond_probs <- sweep(cond_probs, MARGIN = c(3), sigmas_new^2, "/")
        mean_est <- apply(norm_cond_probs, MARGIN = c(1,2),
                          function(p_vec){
                            sum(p_vec*mus_new)
                          })

        mean_ys <- sweep(norm_cond_probs, MARGIN = c(1,2), Y, "*")
        mean_ys <- apply(mean_ys, MARGIN = c(1,2), sum)
        psi_inv <-  apply(norm_cond_probs, MARGIN = c(1,2), sum)
        pred <- qr.fitted(q, as.vector((mean_ys - mean_est)/psi_inv)[!is.na(Y)])
        S <- Y
        S[!is.na(Y)] <- pred
      }
    }
    e <- Y - S

    # Check differences
    dist <- max(c(abs(mus_old - mus_new),
                  abs(sigmas_old - sigmas_new)))

    iter <- iter + 1
    if(verbose) {cat("\r EM iteration:", iter)}
    mus_old <- mus_new
    sigmas_old <- sigmas_new
  }

  if(verbose) {cat("\r Finished with ",iter, "iterations. \n")}

  df_par <- data.frame(mu = mus_old, sigma = sigmas_old)
  rownames(df_par) = 0:C

  return(list(par = df_par,
              fixed = S,
              Z_pred = Z,
              predicted = apply(Z, c(1,2), function(z) mus_old[z+1]) + S,
              iterations = iter))
}
