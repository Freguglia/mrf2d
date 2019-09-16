#' @name dplot
#' @author Victor Freguglia
#' @title Plotting functions for lattice data
#'
#' @description `dplot()` and `cplot()` are functions for plotting lattice data.
#' They are an alternative to base R's `image()` function using `ggplot2`
#' instead.
#' `dplot` is used for discrete data and `cplot` for continuous data, they only
#' differ in the fact that pixel values are treated as a factor in `dplot`,
#' therefore, a discrete scale is used.
#'
#'
#' @param Z A `matrix` object with integers only.
#' @param Y A `matrix` object with continuous values.
#' @param legend `logical` indicating whether a legend should be included or not.
#'
#' @return a `ggplot` object.
#
#' @details Since returns a `ggplot` object, other layers can be added to it
#' using the usual `ggplot2` syntax in order to modify any aspect of the plot.
#'
#' The data frame used to create the object has columns named `x`, `y` and
#' `value`, which are mapped to `x`, `y` and `fill`, respectively, used
#' with `geom_tile()`.
#'
#' @examples
#' # Plotting discrete data
#' dplot(Z_potts)
#'
#' #Making it continuous
#' cplot(Z_potts + rnorm(length(Z_potts)))
#'
#' #Adding extra ggplot layers
#' library(ggplot2)
#' dplot(Z_potts) + ggtitle("This is a title")
#' dplot(Z_potts, legend = TRUE) + scale_fill_brewer(palette = "Set1")
#'
#'
#' @import ggplot2
#' @export
dplot <- function(Z, legend = FALSE){
  df <- data.frame(x = as.vector(row(Z)),
                   y = as.vector(col(Z)),
                   value = as.vector(Z))
  p <- ggplot(df, aes_string(x = "x", y = "y", fill = "factor(value)")) +
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_grey(start = 0, end = 0.9, na.value = "white") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(margin = margin(l=5)))
  if(!legend) { p <- p + theme(legend.position = "none")}
  return(p)
}

#' @rdname dplot
#'
#' @export
cplot <- function(Y, legend = TRUE){
  df <- data.frame(x = as.vector(row(Y)),
                   y = as.vector(col(Y)),
                   value = as.vector(Y))
  p <- ggplot(df, aes_string(x = "x", y = "y", fill = "value")) +
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient(low = "black", high = "gray90", na.value = "white") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(margin = margin(l=5)))
  if(!legend) { p <- p + theme(legend.position = "none")}
  return(p)
}
