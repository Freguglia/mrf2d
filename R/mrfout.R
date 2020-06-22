#' @name mrfout
#' @title MRF fitting functions output
#' @param x a `mrfout` object.
#' @param ... other arguments not used by this method

#' @rdname mrfout
#' @export
print.mrfout <- function(x, ...){
  M <- x$mrfi@Rmat
  pos <- split(M, rep(1:nrow(M), times = 2))
  pos <- lapply(pos, function(x){
    paste0("(", x[1], ",", x[2], ")", collapse = "")
  })
  pos <- unlist(pos)
  cat("Model fitted via", x$method, "\n")
  cat(nrow(M), "interacting positions:", pos, "\n")
  cat("family:", x$family, "\n")
}
