# Converts an array of pairwise counts to a vector of sufficient statistics.
suf_stat <- function(arr, family){
  n_R <- dim(arr)[3]
  C <- dim(arr)[1] - 1
  if(C == 0) stop("Each slice needs at least two rows/columns.")
  if(any(arr < 1 & arr > 0)) stop("'arr' must contain counts instead of proportions.")

  if(family == "onepar"){
    sum((1 - array(diag(C+1), dim = c(C+1, C+1, n_R)))*arr)

  } else if(family == "oneeach"){
    apply(arr, MARGIN = 3, function(m){
      sum((1 - diag(C+1))*m)
    })

  } else if(family == "absdif"){
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(1:C), function(v) {
        sum(m[abs(col(m) - row(m)) == v])
      })
    }))

  } else if(family == "dif"){
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(-C:-1,1:C), function(v) {
        sum(m[col(m) - row(m) == v])
      })
    }))

  } else if(family == "free"){
    as.vector(apply(arr, MARGIN = 3, function(m){
      return(as.vector(m)[-1])
    }))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}
