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

#' @name bold5000
#' @title BOLD5000 neuroimaging data
#'
#' @description An image extracted from the "BOLD5000" open dataset. It was read from
#' the file in path `BOLD5000/DS001499/SUB-CSI2/SES-16/ANAT/SUB-CSI2_SES-16_T1W.NII.GZ`,
#' available at the OpenNeuro platform (https://openneuro.org/datasets/ds001499/versions/1.3.0).
#'
#' @details The file was read using the `oro.nifti` package and the image was extracted from the
#' matrix in slice 160.
#'
#' @references Chang, N., Pyles, J. A., Marcus, A., Gupta, A., Tarr, M. J., & Aminoff, E. M. (2019).
#'  BOLD5000, a public fMRI dataset while viewing 5000 visual images. Scientific data, 6(1), 1-18.
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \url{https://arxiv.org/abs/2006.00383}
#'
#' @docType data
"bold5000"
