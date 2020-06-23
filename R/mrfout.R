#' @name mrfout
#' @title MRF fitting functions output
#' @param x a `mrfout` object.
#' @param object a `mrfout` object.
#' @param ... other arguments not used by this method.

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

#' @rdname mrfout
#' @export
summary.mrfout <- function(object, ...){
  M <- object$mrfi@Rmat
  theta <- object$theta
  family <- object$family
  C <- nrow(theta) - 1
  potentials <- lapply(1:nrow(object$mrfi@Rmat),
                       function(x) smr_array(object$theta[,,x, drop = FALSE], family))
  contributions <- sapply(1:nrow(object$mrfi@Rmat),
                          function(x)
                            sum(abs(cohist(object$Z, object$mrfi)[,,x]*object$theta[,,x]))
                          )
  contributions <- unlist(contributions)
  contributions <- contributions/max(contributions)
  codes <- cut(contributions, breaks = c(0, .25, .5, .75, 1), labels = c(".", "*", "**", "***"))
  pos <- split(M, rep(1:nrow(M), times = 2))
  pos <- lapply(pos, function(x){
    paste0("(", x[1], ",", x[2], ")", collapse = "")
  })
  pos <- unlist(pos)

  if(family != "onepar"){
    if(family == "oneeach"){
      cat(sprintf("%10s", ""), "Interaction for:", "\n")
      cat("Position  ", "Different")
      cat(sprintf("%10s", ""),"Rel. Contribution")
      cat("\n")
    } else if (family == "absdif"){
      cat(sprintf("%10s", ""), "Interaction for abs. difference:", "\n")
      cat("Position  ", sprintf(" %-5s",1:C))
      cat(sprintf("%10s", ""),"Rel. Contribution")
      cat("\n")
    } else if (family == "dif"){
      cat(sprintf("%10s", ""), "Interaction for difference:", "\n")
      cat("Position  ", sprintf(" %-6s",c(-C:-1, 1:C)))
      cat(sprintf("%10s", ""),"Rel. Contribution")
      cat("\n")
    } else if (family == "free"){
      cat(sprintf("%10s", ""), "Interaction for pair:", "\n")
      cat("Position  ", sprintf(" %-6s",paste0("(", rep(0:C, C+1), ",", rep(0:C, each = C+1), ")")[-1]))
      cat(sprintf("%10s", ""),"Rel. Contribution")
      cat("\n")
    }
    for(i in seq_along(pos)){
      cat("",sprintf("%-10s", pos[i]))
      cat(sprintf("%7s", sprintf("%.3f", potentials[[i]])))
      cat(sprintf("%10s",""),
          sprintf("%7s", sprintf("%.3f", contributions[i])),
          as.character(codes[i]))
      cat("\n")
    }
  }
}
