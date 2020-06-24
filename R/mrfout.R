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
  mrfi <- object$mrfi
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

  cat("Model adjusted via", object$method, "\n")
  cat("Image dimension:", dim(object$Z), "\n")
  cat(C+1, "colors, distributed as:\n")
  cat(sprintf("%6s", c(0:C, ifelse(any(is.na(object$Z)), "<NA>", ""))), "\n")
  cat(sprintf("%6s", unname(table(object$Z, useNA = "ifany"))), "\n")
  cat("\n")

  if(object$method == "Pseudolikelihood"){
  }


  if(family == "onepar"){

    equals <- apply(cohist(object$Z, mrfi), function(x) sum(diag(x)), MARGIN = 3)
    diffs <- apply(cohist(object$Z, mrfi), function(x) sum(x) - sum(diag(x)), MARGIN = 3)

    cat("Interaction parameter for different-valued pairs:\n", potentials[[1]], "\n\n")
    cat("Position|", sprintf("%7s", "Equal") ,"Different", "\n")
    for(i in seq_along(pos)){
      cat("", sprintf("%8s", paste0(pos[i], "|", collapse = "")))
      cat(sprintf("%8s", sprintf("%7i", equals[i])))
      cat(sprintf("%8s", sprintf("%7i", diffs[i])))
      cat("\n")
    }
    cat("", sprintf("%8s", "total|"))
    cat(sprintf("%8s", sprintf("%7i", sum(equals))))
    cat(sprintf("%8s", sprintf("%7i", sum(diffs))))
  }
  else {
    if(family == "oneeach"){
      cat("Interactions for different-valued pairs:", "\n")
      cat("Position|", sprintf("%6s", "Value"), " Rel. Contribution")
      cat("\n")
    } else if (family == "absdif"){
      cat("Interactions for abs. differences:", "\n")
      cat("Position|", sprintf("%6s",1:C), " Rel. Contribution")
      cat("\n")
    } else if (family == "dif"){
      cat("Interactions for differences:", "\n")
      cat("Position|", sprintf("%6s",c(-C:-1, 1:C)), " Rel. Contribution")
      cat("\n")
    } else if (family == "free"){
      cat("Interaction for pairs of values:", "\n")
      cat("Position|",
          sprintf("%-6s",paste0("(", rep(0:C, C+1), ",", rep(0:C, each = C+1), ")")[-1]),
          " Rel. Contribution")
      cat("\n")
    }
    for(i in seq_along(pos)){
      cat("", sprintf("%8s", paste0(pos[i],"|", collapse = "")))
      cat("", sprintf("%6s", sprintf("%.3f", potentials[[i]])))
      cat("", sprintf("%6s", sprintf("%.3f", contributions[i])),
          as.character(codes[i]))
      cat("\n")
    }
  }
}

#' @rdname mrfout
#' @export
plot.mrfout <- function(x, ...){
  vecs <- lapply(seq_len(nrow(x$mrfi@Rmat)),
                 function(i) smr_array(x$theta[,,i, drop = FALSE], family = x$family))
  df <- do.call(rbind, vecs)
  C <- nrow(x$theta) - 1
  if(x$family == "oneeach" || x$family == "onepar"){
    colnames(df) <- "Different"
  } else if(x$family == "absdif"){
    colnames(df) <- 1:C
  } else if(x$family == "dif"){
    colnames(df) <- (-C:C)[-(C+1)]
  } else if(x$family == "free"){
    colnames(df) <- paste0("(", rep(0:C, C+1), ",", rep(0:C, each = C+1), ")")[-1]
  }
  df <- cbind(x$mrfi@Rmat, df)
  colnames(df)[1:2] <- c("rx", "ry")
  df <- tidyr::pivot_longer(as.data.frame(df), cols = -c("rx", "ry"))
  max_norm <- max(5, max(df$rx), max(df$ry)) + 0.5
  ggplot(df, aes_string(x = "rx", y = "ry", fill = "value")) +
    geom_tile(color = "black") +
    geom_tile(aes(x = 0, y = 0), fill = "black") +
    geom_tile(aes_string(x = "-rx", y = "-ry"), color = "gray35", fill = "gray95", linetype = "dashed") +
    facet_wrap(~name) +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white") +
    theme_minimal() +
    lims(x = c(-max_norm, max_norm), y = c(-max_norm, max_norm)) +
    ggtitle(x$family) + labs(fill = "theta")
}

