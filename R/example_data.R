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

#' @name data_examples
#' @title Example Data
#' @author Victor Freguglia
#'
#' @description `mrf2d` contains a set of simulated fields to illustrate its
#' usage.
#' \describe{
#'  \item{field1}{A binary field sampled from a sparse interaction structure:
#'  `mrfi(1) + c(4,4)`}
#'  \item{hfield1}{A continuous valued field, obtained by Gaussian mixture driven
#'  by `field1`.}
#' }
#' @docType data
"field1"

#' @rdname data_examples
"hfield1"
