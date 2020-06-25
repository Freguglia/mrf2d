#' @name hmrfout
#' @title MRF fitting functions output
#' @param x a `hmrfout` object.
#' @param object a `hmrfout` object.
#' @param ... other arguments not used by this method.

#' @rdname hmrfout
#' @export
print.hmrfout <- function(x, ...){
  C <- nrow(x$theta) - 1
  cat("Gaussian mixture with", C, "components driven by hidden MRF\n")
  cat("Fitted in", x$iterations, "EM algorithm iterations.\n")
}

#' @rdname hmrfout
#' @export
summary.hmrfout <- function(object, ...){
  C <- object$C
  M <- object$mrfi@Rmat
  pos <- split(M, rep(1:nrow(M), times = 2))
  pos <- lapply(pos, function(x){
    paste0("(", x[1], ",", x[2], ")", collapse = "")
  })
  pos <- unlist(pos)
  cat("Gaussian mixture model driven by Hidden MRF fitted by EM-algorithm.\n")
  cat("Image dimensions:", dim(object$Y), "\n")

  cat("Predicted mixture component table:\n")
  cat(sprintf("%6s", c(0:C, ifelse(any(is.na(object$Z)), "<NA>", ""))), "\n")
  cat(sprintf("%6s", unname(table(object$Z_pred, useNA = "ifany"))), "\n")

  cat("Number of covariates (or basis functions):", length(object$covs), "\n")
  cat("Interaction structure considered:", pos, "\n")

  cat("\n")
  cat("Mixture parameters:\n")
  cat(sprintf("%10s", "Component"), sprintf("%6s", 0:C), "\n")
  cat(sprintf("%10s", "mu"), sprintf("%6s", sprintf("%.2f", object$par$mu)), "\n")
  cat(sprintf("%10s", "sigma"), sprintf("%6s", sprintf("%.2f", object$par$sigma)), "\n")

  cat("\nModel fitted in", object$iterations, "iterations.")
}

#' @rdname hmrfout
#' @export
plot.hmrfout <- function(x, ...){
  plot(dplot(x$Z_pred) + ggtitle("Predicted component"))
  par(ask = TRUE)
  if(any(x$fixed != 0)){
    plot(cplot(x$fixed) + ggtitle("Estimated fixed effect"))
  }
  plot(cplot(x$predicted) + ggtitle("Predicted pixel mean"))
  par(ask = FALSE)
}
