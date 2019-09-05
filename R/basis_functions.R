#' @name basis_functions
#' @author Victor Freguglia
#' @title Creation of basis functions
#'
#' @description `fourier_2d()` and `polynomial_2d()` creates a `list` of basis
#' functions to be used as the fixed effect in \code{\link[=fit_ghm]{fit_ghm}}.
#'
#' @param lattice_dim A length 2 numeric vector with lattice dimensions (N,M)
#'  to be used.
#' @param max_freqs A length 2 numeric vector with maximum frequencies considered
#' (x-axis and y-axis direction, respectively).
#' @param poly_deg A length 2 numeric vector with degrees of the bivariate
#'  polynomial to be considered.
#'
#' @return A `list` of functions.
#'
#' @details `fourier_2d()` is for 2-dimensional Fourier transform.
#'
#' @examples
#' \dontrun{
#' fourier_2d(c(10,10), dim(Z_potts))
#' polynomial_2d(c(3,3), dim(Z_potts))
#' }
NULL

#' @rdname basis_functions
#' @export
fourier_2d <- function(max_freqs, lattice_dim){

  if(any(max_freqs < 0) || length(max_freqs) != 2)  {
    stop("'max_freq' must be length 2 vector of positive integers.
         Read the function documentation for details.")}
  if(any(lattice_dim < 0) || length(lattice_dim) != 2)  {
    stop("'lattice_dim' must be length 2 vector of positive integers.
         Read the function documentation for details.")}

  N <- lattice_dim[1]
  M <- lattice_dim[2]
  n <- 0:max_freqs[1]
  m <- 0:max_freqs[2]
  cbn <- expand.grid(n,m)

  l <- split(cbn, 1:nrow(cbn))[-1]
  l_sin_sin <- lapply(l, function(nm)
    return(
      function(x, y, N, M) sin(2*(x-1)*pi*nm[1,1]/N)*sin(2*(y-1)*pi*nm[1,2]/M)))
  l_cos_sin <- lapply(l, function(nm)
    return(
      function(x, y, N, M) cos(2*(x-1)*pi*nm[1,1]/N)*sin(2*(y-1)*pi*nm[1,2]/M)))
  l_sin_cos <- lapply(l, function(nm)
    return(
      function(x, y, N, M) sin(2*(x-1)*pi*nm[1,1]/N)*cos(2*(y-1)*pi*nm[1,2]/M)))
  l_cos_cos <- lapply(l, function(nm)
    return(
      function(x, y, N, M) cos(2*(x-1)*pi*nm[1,1]/N)*cos(2*(y-1)*pi*nm[1,2]/M)))

  return(c(l_sin_sin, l_cos_cos, l_cos_sin, l_sin_cos))
}

#' @rdname basis_functions
#' @export
polynomial_2d <- function(poly_deg, lattice_dim){
  if(any(poly_deg < 0) || length(poly_deg) != 2)  {
    stop("'max_freq' must be length 2 vector of positive integers.
         Read the function documentation for details.")}
  if(any(lattice_dim < 0) || length(lattice_dim) != 2)  {
    stop("'lattice_dim' must be length 2 vector of positive integers.
         Read the function documentation for details.")}

  N <- lattice_dim[1]
  M <- lattice_dim[2]
  n <- 0:poly_deg[1]
  m <- 0:poly_deg[2]
  cbn <- expand.grid(n,m)

  l <- split(cbn, 1:nrow(cbn))[-1]
  l_poly <- lapply(l, function(nm) return(
    function(x, y, N, M) return((x - N/2)^nm[1,1]*(y - M/2)^nm[1,2]) ))

  return(l_poly)
}

globalVariables(c(".x_basis", ".y_basis"))

# Creates a data frame based on a list of basis functions.
# Simplifies the use in other functions.
#' @importFrom stats sd
basis_function_df <- function(fn_list, N, M, standardize = TRUE){
  df <- reshape2::melt(matrix(0, N, M), varnames = c(".x_basis",".y_basis"))[ ,1:2]
  for(i in seq_along(fn_list)){
    df <- dplyr::mutate(df, fn_list[[i]](.x_basis, .y_basis, N, M))
    colnames(df)[i+2] <- paste0("f",i)
  }

  if(standardize){
    for(i in seq_along(fn_list)){
      df[ ,(i+2)] <- (df[ ,(i+2)] - mean(df[ ,(i+2)]))/sd(df[ ,(i+2)])
    }
  }

  df <- df[,!is.na(df[1,])]
  return(df)
}
