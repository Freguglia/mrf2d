# Converts an array of pairwise counts to a vector of sufficient statistics.
suf_stat <- function(arr, family){
  n_R <- dim(arr)[3]
  C <- dim(arr)[1] - 1
  if(C == 0) stop("Each slice needs at least two rows/columns.")
  if(any(arr < 1 & arr > 0)) stop("'arr' must contain counts instead of proportions.")

  if(family == "onepar"){
    sum((1 - array(diag(C+1), dim = c(C+1, C+1, n_R)))*arr)

  } else if(family == "oneeach"){
    apply(arr, MARGIN = 3, function(m){
      sum((1 - diag(C+1))*m)
    })

  } else if(family == "absdif"){
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(1:C), function(v) {
        sum(m[abs(col(m) - row(m)) == v])
      })
    }))

  } else if(family == "dif"){
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(-C:-1,1:C), function(v) {
        sum(m[col(m) - row(m) == v])
      })
    }))

  } else if(family == "free"){
    as.vector(apply(arr, MARGIN = 3, function(m){
      return(as.vector(m)[-1])
    }))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

#' @name smr_stat
#' @author Victor Freguglia
#' @title Summary Statistics
#'
#' @description Computes the summary count statistics of a field given an
#' interaction structure and a restriction family.
#'
#' @details The order the summarized counts appear in the summary vector matches
#' the order in \code{\link[=smr_array]{smr_array()}}.
#'
#' @inheritParams fit_pl
#'
#' @return A numeric vector with the summarized counts.
#'
#' @examples
#' smr_stat(Z_potts, mrfi(1), "onepar")
#' smr_stat(Z_potts, mrfi(1), "oneeach")
#'
#' @export
smr_stat <- function(Z, mrfi, family){
  C <- max(Z)
  smr_array <- table_relative_3d(Z, mrfi@Rmat, C)
  return(suf_stat(smr_array, family))
}

#' @name smr_array
#' @author Victor Freguglia
#' @title Summarized representation of theta arrays
#'
#' @description Creates a vector with only the free parameters from an array.
#'
#' @inheritParams fit_pl
#' @inheritParams rmrf2d
#'
#' @details The order the parameters appear in the vector matches
#' the order in \code{\link[=smr_stat]{smr_stat()}}.
#'
#' @return A numeric vector with the free parameters of `theta`.
#'
#' @examples
#' smr_array(theta_potts, "onepar")
#' smr_array(theta_potts, "oneeach")
#'
#' @export
smr_array <- function(theta, family){
  array_to_vec(theta, family)
}
