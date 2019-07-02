#' @rdname mrfi-class
#' @title The `mrfi` class
#' @description  A representation of the interaction structure for a spatially-stationary
#' Markov Random Field.
#'
#' @slot Rmat A 2-column `matrix` where each row represents a relative position of
#' interaction.
#' @slot n_neis The number of intracting positions.
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
#' @param max_norm a `numeric` value. All points with norm \eqn{\le} `max_dist` are included.
#' @param norm_type a `character` indicating the norm type used. Possible values are "m", "1", "2", etc. See \code{\link[=norm]{norm}} for details.
#' @param positions a `list` of `numeric` vectors of length 2. Each vector corresponds to a relative position included.
#' @return A \code{\link[=mrfi-class]{mrfi}} object.
#' @importFrom methods new
#' @export
mrfi <- function(max_norm = 1, positions = NULL, norm_type = "1"){
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
  new("mrfi", Rmat = df, n_neis = nrow(df))
}


#' @rdname plot.mrfi
#' @aliases mrfi-plot
#' @title Plotting of `mrfi` objects.
#'
#' @description Plots a visual description of the interaction structure
#'  described in a `mrfi` object. The black tile represents a reference pixel
#'  and gray tiles are shown in relative positions where pixels are
#'  conditionally dependent.
#'
#' @param x A \code{\link[=mrfi-class]{mrfi}} object.
#' @param no_axis `logical` value indicating whether the axis and grid lines are present.
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
