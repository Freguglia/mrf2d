#' @name Z_potts
#' @title Example objects from `mrf2d`
#' @author Victor Freguglia
#'
#' @description `Z_potts` and `theta_potts` are example objects for `mrf2d`.
#'
#' `Z_potts` is a `matrix` object containing an observed lattice of a 3 color
#' (C = 2) Potts model.
#'
#' `theta_potts` is the parameter array used to sample it,
#' it consists of a configuration with one parameter (-1.0) and two relative
#' positions (to be used with a nearest-neighbor structure).
#'
#' @examples
#' theta_potts
#' dplot(Z_potts)
#'
#' @docType data
NULL

#' @name theta_potts
#' @rdname Z_potts
NULL
