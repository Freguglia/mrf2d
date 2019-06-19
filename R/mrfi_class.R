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
            cat(" rx  | ry  \n")
            for(i in seq_len(min(5,object@n_neis))){
              cat(" ",object@Rmat[i,1],"  | ", object@Rmat[i,2], "\n")
            }
            if(object@n_neis > 5) {
              cat("... and", (object@n_neis-5), "more.")
            }
          })

#' @name plot.mrfi
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
#' @importFrom ggplot2 ggplot aes_string geom_tile theme_minimal theme_void lims
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

            max_dist <- max(5, max(df$rx), max(df$ry)) + 0.5
            p <- ggplot(df, aes_string(x = "rx", y = "ry")) +
              geom_tile(fill = "gray", color = "black") +
              geom_tile(data = df_center, fill = "black") +
              theme_minimal() +
              lims(x = c(-max_dist, max_dist), y = c(-max_dist, max_dist))
            if(no_axis) {p <- p + theme_void()}
            p
          })
