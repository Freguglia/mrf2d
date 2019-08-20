#' @name mrfi-class
#' @title An S4 class to represent interactions for MRFs in 2d lattices
#' @description  The `mrfi` S4 class is a representation of the interaction
#'  structure for a spatially-stationary Markov Random Field.
#'
#' The function \code{\link[=mrfi]{mrfi()}} provides an interface for creation
#'  `mrfi` objects. A `plot` method is also available for visualization.
#'
#' @slot Rmat A 2-column `matrix` where each row represents a relative position
#' of interaction.
#' @slot n_neis The number of intracting positions.
#'
#' @details The interaction structure is defined by the list of relative
#' positions in it. For a specific pixel, conditional to the values of pixels in
#' these relative positions from it, its value is independent of any other pixel
#' in the image.
#'
#' The relative positions are indentified by two integers `rx` and `ry`
#' representing the "shift" in the `x`-axis and `y`-axis respectively. As an
#' example: The relative position `(1,0)` representes the pixel in the immediate
#' right position, while `(-1,0)` the left one.
#'
#' Note that the inclusion of a relative position to the dependence also implies
#' its opposite direction is not conditionally independent (commutativeness of
#' dependence), but only one is actually included to the `mrfi` object.
#'
#' To illustrate that, a nearest neighbor dependence structure can be specified
#' by:
#'
#' \code{mrfi(1)}
#'
#' Note that it only includes the positions `(1,0)` and `(0,1)`, but when
#' visualizing it, for example, `mrf2d` understands the opposite directions
#' are also conditionally dependent, as in
#'
#' \code{plot(mrfi(1))}.
#'
#' @examples
#' plot(mrfi(max_norm = 2, norm_type = "1"))
#' plot(mrfi(max_norm = 2, norm_type = "m"))
#' plot(mrfi(max_norm = 2, norm_type = "1", positions = list(c(4,4))))
#'
#' @exportClass mrfi
setClass("mrfi",
         representation(Rmat = "matrix",
                        n_neis = "numeric"))

setMethod("show", "mrfi",
          function(object){
            cat(object@n_neis, "interacting positions.\n")
            cat("  rx     ry")
            for(i in seq_len(min(5,object@n_neis))){
              cat( "\n"," ",object@Rmat[i,1],"    ", object@Rmat[i,2])
            }
            if(object@n_neis > 5) {
              cat("  ... and", (object@n_neis-5), "more.")
            }
          })


#' @name mrfi
#' @title Creation of \code{\link[=mrfi-class]{mrfi}} objects.
#' @author Victor Freguglia
#'
#' @description `mrfi()` creates an object of class `mrfi` based on a distance
#' rule and optionally a list of relative positions. The argument `max_norm` and
#' `norm_type` can be used to automatically include all positions within a
#' "range" defined by the norm type chosen and distance using that norm.
#'
#' A list of relative positions may also be included to specify sparse
#' interaction structures, for example.
#'
#' @param max_norm a `numeric` value. All points with norm \eqn{\le} `max_dist`
#'  are included.
#' @param norm_type a `character` indicating the norm type used. Possible values
#'  are "m", "1", "2", etc. See \code{\link[=norm]{norm}} for details.
#' @param positions a `list` of `numeric` vectors of length 2. Each vector
#'  corresponds to a relative position included.
#'
#' @return A \code{\link[=mrfi-class]{mrfi}} object.
#'
#' @note If a position in `positions` is already included due to the
#' `max_norm` and `norm_type` specification, the second ocurrence is ignored.
#' The same is valid for repeated or opposite positions in `positions`.
#'
#' @examples
#' mrfi(1)
#' mrfi(2)
#' mrfi(2, norm_type = "m")
#' mrfi(1, positions = list(c(4,4), c(-4,4)))
#'
#' #Repeated positions are handled automatically
#' mrfi(1, positions = list(c(1,0), c(2,0)))
#'
#' @importFrom methods new
#' @export
mrfi <- function(max_norm = 1, norm_type = "1", positions = NULL){
  if(max_norm < 0){stop("'max_norm' must be greater than or equal 0.")}
  if(max_norm > 0){
    df <- expand.grid(x = -max_norm:max_norm, y = 0:max_norm)
    df <- df[apply(df, MARGIN = 1,
                   function(m) norm(matrix(m), type = norm_type)) <= max_norm,]
  } else {
    df <- data.frame(x = c(0,0), y = c(0,0))
  }

  if(!is.null(positions)){
    if(!is.list(positions)){
      stop("'positions' must be a list of relative positions.")
    } else if(any(!unlist(lapply(positions, is.numeric)))) {
      stop("'positions' must be a list of relative positions (numeric).")
    } else {
      df <- rbind(as.matrix(df), do.call(rbind, positions))
    }
  }

  df_minus <- -df
  str_vec_df <- apply(df, MARGIN = 1, paste0, collapse = "__")
  str_vec_minus <- apply(df_minus, MARGIN = 1, paste0, collapse = "__")
  to_remove <- which(str_vec_df %in% str_vec_minus)


  while(length(to_remove) > 0){
    if(length(to_remove) > 1) {
      df <- df[-to_remove[1], ]
      df_minus <- -df
      str_vec_df <- apply(df, MARGIN = 1, paste0, collapse = "__")
      str_vec_minus <- apply(df_minus, MARGIN = 1, paste0, collapse = "__")
    } else {
      df <- df[-to_remove,]
      break
    }
    to_remove <- which(str_vec_df %in% str_vec_minus)
  }
  df <- matrix(unlist(df), ncol = 2)
  df <- unique(df)
  new("mrfi", Rmat = df, n_neis = nrow(df))
}


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
#' @param no_axis `logical` value indicating whether the axis and grid lines
#'  are used. If `TRUE` it simply adds `theme_void()` to the `ggplot` object.
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
#' @exportMethod plot
setMethod("plot", signature(x = "mrfi", y = "missing"),
          definition = function(x, no_axis = FALSE){
            df <- as.data.frame(x@Rmat)
            names(df) <- c("rx", "ry")
            df2 <- df
            df2$rx <- -df$rx
            df2$ry <- -df$ry
            df <- rbind(df,df2)

            df_center <- data.frame(rx = 0, ry = 0)

            max_norm <- max(5, max(df$rx), max(df$ry)) + 0.5
            p <- ggplot(df, aes_string(x = "rx", y = "ry")) +
              geom_tile(fill = "gray", color = "black") +
              geom_tile(data = df_center, fill = "black") +
              theme_minimal() +
              lims(x = c(-max_norm, max_norm), y = c(-max_norm, max_norm))
            if(no_axis) {p <- p + theme_void()}
            p
          })
