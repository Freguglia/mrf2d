## These functions convert arrays of potentials to vector of free parameters.
## Used in maximum pseudo-likelihood estimation, where optim() requires a vector.

# Checks if an array is a valid potential definition.
is_valid_array <- function(arr, family){
  n_r <- dim(arr)[3]
  C <- dim(arr)[1] - 1
  if(dim(arr)[2] != (C+1)) {
    stop("Invalid array: number of rows and columns differ")}
  if(C < 1) stop("At least two different possible values are required.")

  if(family == "onepar"){
    b <- arr[2,1,1]
    diags_ok <- all(apply(arr, MARGIN = 3, function(m) all(diag(m) == 0)))
    rest_ok <- all(apply(arr, MARGIN = 3, function(m) all(m[row(m)!=col(m)] == b)))
    return(diags_ok & rest_ok)

  } else if(family == "oneeach") {
    all(apply(arr, MARGIN = 3, function(m) {
      b <- m[1,2]
      diags_ok <- all(diag(m) == 0)
      rest_ok <- all(m[row(m)!=col(m)] == b)
      return(diags_ok & rest_ok)
    }))

  } else if(family == "absdif") {
    all(apply(arr, MARGIN = 3, function(m) {
      one_by_dif <- sapply(0:C, function(v){
        length(unique(m[abs(col(m)-row(m))==v])) == 1
      })
      zero_diags <- all(apply(arr, MARGIN = 3, function(m) all(diag(m) == 0)))
      return(all(one_by_dif) & zero_diags)
    }))

  } else if(family == "dif") {
    all(apply(arr, MARGIN = 3, function(m){
      #Check if matrix is valid:
      # Are potentials equal for each difference?
      one_by_dif <- sapply(-C:C, function(v){
        length(unique(m[col(m)-row(m)==v])) == 1
      })
      # are diagonals zero?
      zero_diags <- all(apply(arr, MARGIN = 3, function(m) all(diag(m) == 0)))
      return(all(one_by_dif) & zero_diags)
    }))

  } else if(family == "free") {
    all(apply(arr, MARGIN = 3, function(m) {
      m[1,1] == 0
    }))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts vector to appropriate array of potentials.
vec_to_array <- function(vec, family, C, n_R){
  if(family == "onepar") {
    if(length(vec) != 1) { stop("'vec' must have length 1 for family 'onepar'.") }
    sanitize_theta(array(1 - diag(C+1), dim = c(C+1, C+1, n_R)) * vec)

  } else if(family == "oneeach"){
    if(length(vec) != n_R) { stop("'vec' must have length equal to 'n_R' for family 'oneeach'.") }
    sanitize_theta(array(1 - diag(C+1), dim = c(C+1, C+1, n_R)) * rep(vec, each = (C+1)^2))

  } else if(family == "absdif") {
    if(length(vec) != n_R*C) { stop("'vec' must have length n_R*C for family 'absdif'.")}
    sanitize_theta(simplify2array(lapply(1:n_R, function(i){
      v <- vec[(1+(C*(i-1))):(C*i)]
      m <- diag(C+1)*0
      for(j in 1:C){
        m[abs(col(m) - row(m)) == j] <- v[j]
      }
      return(m)
    })))

  } else if(family == "dif") {
    # Potential associated with zero difference (diagonal) is the sum of others.
    # The order is from -C to -1 then 1 to C.
    if(length(vec) != n_R*2*C) { stop("'vec' must have length n_R*2*C for family 'dif'.")}
    sanitize_theta(simplify2array(lapply(1:n_R, function(i){
      v <- vec[(1+(2*C*(i-1))):(2*C*i)]
      m <- diag(C+1)*0
      k <- 1
      for(j in c(-C:-1,1:C)){
        m[col(m) - row(m) == j] <- v[k]
        k <- k+1
      }
      return(m)
    })))

  } else if(family == "free") {
    if(length(vec) != (n_R*((C+1)^2 - 1))) { stop("'vec' must have length n_R*((C+1)^2 - 1) for family 'free'.")}
    # Parameters are filled by Columns.
    sanitize_theta(simplify2array(lapply(1:n_R, function(i) {
      matrix(c(0, vec[(1 + (i-1)*(((C+1)^2)-1)) : (i*((C+1)^2) - i)]),
             nrow = C+1)
    })))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts array of potentials to a vector
array_to_vec <- function(arr, family){
  if(!is_valid_array(arr, family)){ stop(paste0("Not a valid array for family '", family, "'"))}
  n_R <- dim(arr)[3]
  C <- dim(arr)[1] - 1
  if(C == 0) stop("Each slice needs at least two rows/columns.")

  if(family == "onepar") {
    return(arr[2,1,1])

  } else if(family == "oneeach") {
    return(apply(arr, MARGIN = 3, function(m) {
      m[1,2]
    }))

  } else if(family == "absdif") {
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(1:C), function(v) {
        m[abs(col(m) - row(m)) == v][1]
      })
    }))

  } else if(family == "dif") {
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(-C:-1,1:C), function(v) {
        m[col(m) - row(m) == v][1]
      })
    }))

  } else if(family == "free") {
    as.vector(apply(arr, MARGIN = 3, function(m){
      return(as.vector(m)[-1])
    }))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts matrix to one 1-slice array (used in case n_R = 1 and user passes
# a matrix by accident)
matrix_to_array <- function(m) { return(array(m, dim = c(dim(m),1))) }
sanitize_theta <- function(theta) {
  if(is.matrix(theta)) {theta <- matrix_to_array(theta)}
  rownames(theta) <- colnames(theta) <- 0:(dim(theta)[1]-1)
  return(theta)
}
