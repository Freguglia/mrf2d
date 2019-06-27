#' @rdname dplot
#' @aliases cplot
#' @author Victor Freguglia
#' @title Elegant image plots based on matrices
#'
#' @description An alternative to base R's `image` function using `ggplot2`.
#' Use `dplot` for discrete valued matrices and `cplot` for continuous value.
#'
#' @param Z A `matrix` object with integers only.
#' @return a `ggplot` object.
#' @import ggplot2
#' @export
dplot <- function(Z, legend = FALSE){
  df <- data.frame(x = as.vector(row(Z)),
                   y = as.vector(col(Z)),
                   value = as.vector(Z))
  ggplot(df, aes_string(x = "x", y = "y", fill = "factor(value)")) +
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_grey(start = 0, end = 0.9) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(margin = margin(l=5)))
}

#' @rdname dplot
#' @name cplot
#' @param Y A `matrix` object with continuous values.
#' @export
cplot <- function(Y, legend = TRUE){
  df <- data.frame(x = as.vector(row(Y)),
                   y = as.vector(col(Y)),
                   value = as.vector(Y))
  ggplot(df, aes_string(x = "x", y = "y", fill = "value")) +
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient(low = "black", high = "gray90") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(margin = margin(l=5)))
}
