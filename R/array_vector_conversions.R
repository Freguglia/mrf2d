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
        length(unique(m[row(m)-col(m)==v])) == 1
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
    return(array(diag(C+1)*2 - 1, dim = c(C+1, C+1, n_R))*vec)
  } else if(family == "dif") {
    if(length(vec) != n_R*2*C) {
      stop("'vec' must have length n_R*2*C for family 'dif'.")}

  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts array of potentials to a vector
array_to_vec <- function(arr, family){
  if(family == "onepar") {
    stopifnot(is_valid_array(arr, family))
    return(arr[1,1,1])
  } else {
    stop("'", family, "' is not an implemented family.")
  }
}

# Converts matrix to one 1-slice array (used in case n_R = 1 and user passes
# a matrix by accident)
matrix_to_array <- function(m) { return(array(m, dim = c(dim(m),1))) }
sanitize_theta <- function(theta) {
  if(is.matrix(theta)) {return(matrix_to_array(theta))}
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

