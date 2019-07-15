#' @name fit_ghm
#' @author Victor Freguglia
#' @title EM for Gaussian Hidden Markov field
#'
#' @description Em algorithm fitting of Gaussian mixture driven by a hidden
#'  Markov random field with known parameters.
#'
#' @details Add details about the estimation procedure.
#'
#' @param Y A matrix of observed pixel values.
#' @inheritParams pl_mrf2d
#' @param fixed_fn A list of functions `fn(x,y)` to be considered as a fixed effect.
#' @param maxiter The maximum number of iterations allowed. Defaults to 100.
#' @param equal_vars `logical` indicating if the mixture model has the same
#'  variances in all components
#' @param icm_cycles Number of steps used in the Iterated Conditional Modes
#'  algorithm executed in each interaction. Defaults to 6.
#'
#' @return A list containing some stuff to be documented.
#' @export
fit_ghm <- function(Y, mrfi, theta, fixed_fn = list(),
                    equal_vars = TRUE, init_mus = NULL,
                    init_sigmas = NULL,
                    maxiter = 100, max_dist = 10^-3,
                    icm_cycles = 6, verbose = TRUE){

  Rmat <- mrfi@Rmat

  # check if theta's dimensions matches R's
  if(dim(theta)[3] != nrow(Rmat)) {stop("Invalid dimensions.")}
  C <- dim(theta)[1] - 1

  # compute initial fixed effect
  N <- nrow(Y); M <- ncol(Y)
  if(length(fixed_fn) > 0){
    df_fixed <- basis_function_df(fixed_fn, N, M, standardize = TRUE)
    fit_reg <- lm(gl ~ 0 + ., data = cbind(df_fixed[,-(1:2)], gl = as.vector(Y - mean(Y))))
    S <- matrix(predict(fit_reg), nrow = N, ncol = M)
    e <- Y - S
  } else {
    e <- Y
    S <- Y*0
  }

  # Initialize variables
  if(is.null(init_mus) | is.null(init_sigmas)) {mus_old <- seq(min(e), max(e), length.out = C+1)
    if(verbose) cat("\r Fitting independent mixture to obtain initial parameters. \n")
    ind_fit <- fit_ghm(e, mrfi, theta*0, fixed_fn, equal_vars,
                       init_mus = seq(min(e), max(e), length.out = C+1),
                       init_sigmas = rep(diff(range(e))/(2*C), C+1),
                       maxiter, max_dist, icm_cycles, verbose = FALSE)
    mus_old <- ind_fit$par$mu
    sigmas_old <- ind_fit$par$sigma
  } else {
    mus_old <- init_mus
    sigmas_old <- init_sigmas
  }
  Z <- matrix(sample(0:C, size = N*M, replace = TRUE), nrow = nrow(Y))

  # Iterate
  iter <- 0; dist <- Inf
  while(iter < maxiter && dist > max_dist){

    # Compute MAP estimates of Z
    Z <- icm_gaussian_cpp(e, Rmat, Z, theta, mus_old, sigmas_old, icm_cycles)

    # Compute conditional probabilities
    cond_probs <- cprob_ghm_all(Z, Rmat, theta, mus_old, sigmas_old, e)

    # Update mus, sigmas and S
    ## mus and sigmas
    if(equal_vars){
      mus_new <- apply(cond_probs, MARGIN = 3, function(x) sum(x*e)/sum(x))
      sigmas_new <- rep(sqrt(sum(simplify2array(lapply(mus_new,
                                                       function(mu) (e - mu)^2))*cond_probs)/(N*M)), C+1)
    }
    else {stop("Different variances still not implemented.")}

    ## update S
    if(length(fixed_fn) > 0){
      fit_reg <- lm(gl ~ 0 + .,
                    data = cbind(df_fixed[,-(1:2)],
                                 gl = as.vector(Y - apply(cond_probs, MARGIN = c(1,2),
                                                          function(p) sum(p*mus_new)))))
      S <- matrix(predict(fit_reg), nrow = N, ncol = M)
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
