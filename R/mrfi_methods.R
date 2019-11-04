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

#' @rdname mrfi-class
#'
#' @param x `mrfi` object.
#'
#' @exportMethod as.list
setMethod("as.list", signature(x = "mrfi"),
          definition = function(x){
            unname(split(x@Rmat, rep(1:nrow(x@Rmat), ncol(x@Rmat))))
          })

#' @rdname mrfi-class
#'
#' @param i vector of indexes to extract interacting positions.
#'
#' @return `as.list()` converts the `mrfi` object to a list of interacting
#' positions (length 2 vectors).
#'
#' `[[` converts to list and subsets it.
#'
#' `[` subsets the `mrfi` object and returns another `mrfi` object.
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
