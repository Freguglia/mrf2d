#' @rdname plot.mrfi
#' @aliases mrfi-plot
#' @title Plotting of `mrfi` objects.
#' @author Victor Freguglia
#'
#' @description Plots a visual representation of the interaction structure
#'  described in a `mrfi` object. The black tile represents a reference pixel
#'  and gray tiles are shown in relative positions with dependent pixels.
#'
#'  A `ggplot` object is used, therefore, the user can load the `ggplot2`
#'  package and add more `ggplot` layers to freely customize the `plot`.
#'
#' @param x A \code{\link[=mrfi-class]{mrfi}} object.
#' @param include_axis `logical` indicating whether the axis and grid lines
#'  are included. If `FALSE` `theme_void()` is added to the `ggplot` object.
#' @param include_opposite ´logical` whether opposite directions should be
#'  included in the visualization of the dependence structure.
#' @param ... other arguments not used by this method.
#'
#' @return A `ggplot` object using `geom_tile()` to represent interacting
#' relative positions.
#'
#' @details The `data.frame` used for the `ggplot` call has columns names `rx`
#' and `ry` repŕesenting the relative positions.
#'
#' @examples
#' plot(mrfi(1))
#'
#' library(ggplot2)
#' plot(mrfi(1)) + geom_tile(fill = "red")
#' plot(mrfi(1)) + geom_tile(fill = "blue") + theme_void()
#'
#' plot(mrfi(1)) + geom_text(aes(label = paste0("(",rx,",",ry,")")))
#'
#' @export
plot.mrfi <- function(x, include_axis = FALSE,
                      include_opposite = TRUE,
                      ...){
  df <- as.data.frame(x@Rmat)
  names(df) <- c("rx", "ry")
  df_center <- data.frame(rx = 0, ry = 0)

  max_norm <- max(5, max(df$rx), max(df$ry)) + 0.5
  p <- ggplot(df, aes_string(x = "rx", y = "ry")) +
    geom_tile(fill = "gray", color = "black") +
    geom_tile(data = df_center, fill = "black") +
    theme_minimal()
  if(include_opposite){p <- p +
    geom_tile(data = data.frame(rx = -df$rx, ry = -df$ry),
              linetype = "dashed", color = "gray55",
              fill = "gray95")}
  if(!include_axis) {p <- p + theme_void()}
  p + lims(x = c(-max_norm, max_norm), y = c(-max_norm, max_norm))
}

#' @rdname mrfi-class
#'
#' @param x `mrfi` object.
#' @param ... other arguments not used by this method.
#'
#' @return `as.list()`: converts the `mrfi` object to a list of interacting
#' positions (list of length-2 vectors).
#'
#' @export
as.list.mrfi <- function(x, ...){
  unname(split(x@Rmat, rep(1:nrow(x@Rmat), ncol(x@Rmat))))
}

#' @rdname mrfi-class
#' @exportMethod length
setMethod("length", signature(x = "mrfi"),
          definition = function(x){
            nrow(x@Rmat)
          })

#' @rdname mrfi-class
#'
#' @param i vector of indexes to extract interacting positions.
#'
#' @return `[[`: converts to list and subsets it.
#'
#' `[`: subsets the `mrfi` object and returns another `mrfi` object.
#'
#' `+`: computes the union of the interaction structure in a `mrfi` object with
#' a `numeric` representing an additional position to include or another `mrfi`
#' object.
setMethod("[[", signature = c("mrfi", "numeric", "missing"),
          definition = function(x, i){
            m <- x@Rmat[i,,drop = FALSE]
            unname(split(m, rep(1:nrow(m), ncol(m))))
          })

#' @rdname mrfi-class
setMethod("[", signature = c("mrfi", "numeric", "missing"),
          definition = function(x, i){
            m <- x@Rmat[i,,drop = FALSE]
            new("mrfi", Rmat = m)
          })

mrfi_union <- function(mrfi1, mrfi2){
  return(mrfi(0, positions = union(as.list(mrfi1), as.list(mrfi2))))
}

mrfi_diff <- function(mrfi1, mrfi2){
  return(mrfi(0, positions =
                setdiff(as.list(mrfi1),
                        c(as.list(mrfi2), lapply(as.list(mrfi2), '*', -1L)))))
}

#' @rdname mrfi-class
#'
#' @description Simple operations are provided to include (set union)
#' new interacting positions to a `mrfi` object with the `'+'` operator and
#' remove positions (set difference) with `-`. Individual positions can be
#' included/excluded using length-2 vectors in the right hand side. Union and
#' set difference of complete structures can also be computed by adding or
#' subtracting two `mrfi` objects.
#'
#' These operations deal with opposite directions filtering to avoid redundancy
#' in the interaction structure.
#'
#' @param e1,mrfi A `mrfi` object.
#' @param e2 Either a second `mrfi` object or a length 2 `numeric` with the new
#' relative position to include (`+`) or remove (`-`).
#'
#' @examples
#' mrfi(1) + c(2,0)
setMethod("+", signature = c("mrfi", "numeric"),
          definition = function(e1, e2){
            if(length(e2) != 2){
              stop("Right hand side must be a length 2 vector representing a relative position.")
            } else if (any(as.integer(e2) != e2)){
              stop("Right hand side must be a vector with two integer values.")
            } else if(any(sapply(as.list(e1), function(pos){
              all(pos == e2) | all(pos == (-e2))}
              )) & length(e1) > 0){
                return(e1)
            }
            result <- mrfi_union(e1, list(e2))
            return(result)
          })

#' @rdname mrfi-class
#'
#' @examples
#' mrfi(1) - c(1,0)
setMethod("-", signature = c("mrfi", "numeric"),
          definition = function(e1, e2){
            if(length(e2) != 2){
              stop("Right hand side must be a length 2 vector representing a relative position.")
            } else if (any(as.integer(e2) != e2)){
              stop("Right hand side must be a vector of two integers.")
            } else {
              if(length(e1) == 0) return(e1)
              e2 <- as.integer(e2)
              return(mrfi_diff(e1, list(e2, -e2)))
            }
          })

#' @rdname mrfi-class
#'
#' @examples
#' mrfi(1) + mrfi(0, positions = list(c(2,0)))
setMethod("+", signature = c("mrfi", "mrfi"),
          definition = function(e1, e2){
            return(mrfi_union(e1,e2))
          })

#' @rdname mrfi-class
#'
#' @examples
#' mrfi(2) - mrfi(1)
setMethod("-", signature = c("mrfi", "mrfi"),
          definition = function(e1, e2){
            if(length(e1) == 0) return(e1)
            return(mrfi_diff(e1,e2))
          })


#' @rdname mrfi-class
#' @export
mrfi_to_string <- function(mrfi){
  s <- mrfi_to_char(mrfi)
  return(paste0("{", paste(s, collapse = ","), "}"))
}
