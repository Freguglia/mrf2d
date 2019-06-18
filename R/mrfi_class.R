#' @rdname mrfi-class
#' @title The `mrfi` class
#' @description  A representation of the interaction structure for a spatially-stationary
#' Markov Random Field.
#'
#' @slot Rmat A 3-column `matrix` where each row represents a relative position of
#' interaction. First and second columns represent the position and third column
#' the interaction type.
#' @slot n_neis The number of intracting positions.
#' @slot n_types The number of types of interaction.
#' @exportClass mrfi
setClass("mrfi",
         representation(Rmat = "matrix",
                        n_neis = "numeric",
                        n_types = "numeric"))

setMethod("show", "mrfi",
          function(object){
            cat(object@n_neis, "interacting positions of", object@n_types, "types: \n")
            cat(" rx  | ry  | type \n")
            for(i in seq_len(min(5,object@n_neis))){
              cat(" ",object@Rmat[i,1],"   ", object@Rmat[i,2], "    ", object@Rmat[i,3], "\n")
            }
            if(object@n_neis > 5) {
              cat("... and", (object@n_neis-5), "more.")
            }
          })

#' @name plot.mrfi
#' @aliases mrfi-plot
#' @title Plotting of `mrfi` objects.
#'
#' @param x A \code{\link[=mrfi-class]{mrfi}} object.
#' @param no_axis `logical` value indicating whether the axis and grid lines are present.
#' @param no_text `logical` value indicating whether interaction types should be included.
#' @importFrom ggplot2 ggplot aes_string geom_tile geom_text theme_minimal theme_void lims
#' @exportMethod plot
setMethod("plot", signature(x = "mrfi", y = "missing"),
          definition = function(x, no_axis = FALSE, no_text = FALSE){
            df <- as.data.frame(x@Rmat)
            names(df) <- c("rx", "ry", "type")
            df2 <- df
            df2$rx <- -df$rx
            df2$ry <- -df$ry
            df <- rbind(df,df2)

            df_center <- data.frame(rx = 0, ry = 0, type = 0)

            max_dist <- max(5, max(df$rx), max(df$ry)) + 0.5
            p <- ggplot(df, aes_string(x = "rx", y = "ry", label = "type")) +
              geom_tile(fill = "gray", color = "black") +
              geom_tile(data = df_center, fill = "black") +
              theme_minimal() +
              lims(x = c(-max_dist, max_dist), y = c(-max_dist, max_dist))
            if(no_axis) {p <- p + theme_void()}
            if(!no_text) {p <- p + geom_text()}
            p
          })
