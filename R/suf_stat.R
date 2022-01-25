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

  } else if(family == "symmetric"){
    as.vector(apply(arr, MARGIN = 3, function(m){
      (t(m)+(m*lower.tri(m)))[lower.tri(m, diag = TRUE)][-1]
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
#'   * `cohist()` computes the co-ocurrence histogram.
#'   * `smr_stat()` computes the co-ocurrence histogram, then converts it into
#'   a vector of sufficient statistics given a \code{\link[=mrf2d-family]{family}} of restrictions.
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
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \doi{10.18637/jss.v101.i08}
#'
#' @export
smr_stat <- function(Z, mrfi, family){
  C <- max(1, max(Z, na.rm = TRUE))
  smr_array <- table_relative_3d(Z, mrfi@Rmat, C)
  return(suf_stat(smr_array, family))
}

#' @rdname smr_stat
#'
#' @return An array representing the co-ocurrence histogram of `Z` in the relative
#' positions contained in `mrfi`. Each row and column corresponds a pair of values
#' in `(0, ..., C)` and each slice corresponds to
#'
#' @examples
#' cohist(Z_potts, mrfi(1))
#'
#' @export
cohist <- function(Z, mrfi){
  C <- max(Z, na.rm = TRUE)
  coh <- table_relative_3d(Z, mrfi@Rmat, C)
  pos_names <- sapply(as.list(mrfi), paste, collapse = ",")
  dimnames(coh) <- list(0:C, 0:C, paste0("(", pos_names, ")"))
  return(coh)
}

#' @name smr_array
#' @author Victor Freguglia
#' @title Summarized representation of theta arrays
#'
#' @description `smr_array` creates a vector containing only the free parameters from an array
#' given a restriction \code{\link[=mrf2d-family]{family}}. `exapand_array` is the reverse
#' operation, expanding a complete array from the vector of sufficient statistics.
#'
#' @inheritParams fit_pl
#' @inheritParams rmrf2d
#'
#' @details The order the parameters appear in the summarized vector matches
#' the order in \code{\link[=smr_stat]{smr_stat()}}.
#'
#' \code{vec_description()} provides a \code{data.frame} object describing
#' which are the relative positions and interactions associated with each
#' element of the summarized vector in the same order.
#'
#' @return `smr_array` returns a numeric vector with the free parameters of `theta`.
#'
#' @examples
#' smr_array(theta_potts, "onepar")
#' smr_array(theta_potts, "oneeach")
#'
#' @seealso
#'
#' A paper with detailed description of the package can be found at
#' \doi{10.18637/jss.v101.i08}
#'
#' @export
smr_array <- function(theta, family){
  array_to_vec(theta, family)
}

#' @rdname smr_array
#'
#' @param theta_vec A `numeric` vector with the free parameters of a potential
#' array. It's dimension depends on the restriction `family`, `C` and the number
#' of interacting positions on `mrfi`.
#' @param C The maximum value of the field.
#'
#' @return `expand_array` returns a three-dimensional `array` of potentials.
#'
#' @examples
#' expand_array(0.99, family = "onepar", mrfi = mrfi(1), C = 2)
#' expand_array(c(0.1, 0.2), family = "oneeach", mrfi = mrfi(1), C = 3)
#'
#' @export
expand_array <- function(theta_vec, family, mrfi, C){
  theta <- vec_to_array(theta_vec, family, C, nrow(mrfi@Rmat))
  pos_names <- sapply(as.list(mrfi), paste, collapse = ",")
  dimnames(theta) <- list(0:C, 0:C, paste0("(", pos_names, ")"))
  return(theta)
}

#' @rdname smr_stat
#'
#' @inheritParams expand_array
#' @return A `data.frame` describing the relative position
#'  and interaction associated with each potential in the vector
#'  form in each row, in the same order.
#'
#' @export
vec_description <- function(mrfi, family, C){
    pos <- apply(mrfi@Rmat, MARGIN = 1, paste0, collapse = ",")
    pos <- paste0("(", pos, ")")
    if(family == "onepar"){
        res <- data.frame(position = as.factor("all"),
                          interaction = as.factor("different"))

    } else if(family == "oneeach"){
        res <- data.frame(position = as.factor(pos),
                          interaction = as.factor("different"))

    } else if(family == "absdif"){
        ints <- paste0("abs.dif. ",1:C)
        res <- data.frame(position = rep(pos, each = C),
                          interaction = rep(ints, times = length(pos)))

    } else if(family == "dif"){
        ints <- paste0("dif. ", c(-C:-1,1:C))
        res <- data.frame(position = rep(pos, each = 2*C),
                          interaction = rep(ints, times = length(pos)))


    } else if(family == "free"){
        arr <- array(dim=c(C+1, C+1, length(pos)))
        ints <- paste0(rep(0:C, times = C+1), ",", rep(0:C, each = C+1))
        ints <- ints[ints != "0,0"]
        res <- data.frame(position = rep(pos, each = (C+1)^2 - 1),
                          interaction = rep(ints, times = length(pos)))

    } else if(family == "symmetric"){
        m <- matrix(nrow = C+1, ncol = C+1)
        ints <- paste0(row(m)[lower.tri(m,TRUE)]-1, ",", col(m)[lower.tri(m, TRUE)]-1)
        ints <- ints[ints != "0,0"]
        res <- data.frame(position = rep(pos, each = (C+1)*(C+2)/2 - 1),
                          interaction = rep(ints, times = length(pos)))


    } else {
        stop("'", family, "' is not an implemented family.")
    }
    return(res)
}

