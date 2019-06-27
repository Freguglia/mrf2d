# Checks if an array is a valid potential definition.
is_valid_array <- function(arr, family){
  n_r <- dim(arr)[3]
  C <- dim(arr)[1] - 1
  if(dim(arr)[2] != (C+1)) {
    stop("Invalid array: number of rows and columns differ")}

  if(family == "onepar"){
    b <- arr[1,1,1]
    diags_ok <- all(apply(arr, MARGIN = 3, function(m) all(diag(m) == b)))
    rest_ok <- all(apply(arr, MARGIN = 3, function(m) all(m[row(m)!=col(m)] == -b)))
    return(diags_ok & rest_ok)

  } else if(family == "dif") {
    all(apply(arr, MARGIN = 3, function(m){
      #Check if matrix is valid:
      # Are potentials equal for each difference?
      one_by_dif <- sapply(-C:C, function(v){
        length(unique(m[col(m)-row(m)==v])) == 1
      })
      # Do values sum zero?
      zero_sum <- abs(sum(m[ ,1]) + sum(m[1, 2:(C+1)])) < 10^-6
      return(all(one_by_dif) & zero_sum)
    }))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts vector to appropriate array of potentials.
vec_to_array <- function(vec, family, C, n_R){
  if(family == "onepar") {
    if(length(vec) != 1) { stop("'vec' must have length 1 for family 'onepar'.") }
    sanitize_theta(array(diag(C+1)*2 - 1, dim = c(C+1, C+1, n_R))*vec)

  } else if(family == "dif") {
    # Potential associated with zero difference (diagonal) is the sum of others.
    # The order is from -C to -1 then 1 to C.
    if(length(vec) != n_R*2*C) {
      stop("'vec' must have length n_R*2*C for family 'dif'.")}
    simplify2array(lapply(1:n_R, function(i){
      v <- vec[(1+(2*C*(i-1))):(2*C*i)]
      m <- -diag(C+1)*sum(v)
      k <- 1
      for(j in c(-C:-1,1:C)){
        m[col(m) - row(m) == j] <- v[k]
        k <- k + 1
      }
      return(m)
    }))

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts array of potentials to a vector
array_to_vec <- function(arr, family){
  stopifnot(is_valid_array(arr, family))
  n_R <- dim(arr)[3]
  C <- dim(arr)[1] - 1

  if(family == "onepar") {
    return(arr[1,1,1])

  } else if(family == "dif") {
    as.vector(apply(arr, MARGIN = 3, function(m) {
      sapply(c(-C:-1,1:C), function(v) {
        m[col(m) - row(m) == v][1]
      })
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

# Creates a vector (with correct length) representing the independent field
zero_vec <- function(mrfi, family){
  n_R <- mrfi@n_neis
  if(family == "onepar") {
    return(0)
  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

