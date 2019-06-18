#' The mrfi class
#'
#' A representation for the interaction structure of a spatially-stationary
#' Markov Random Field.
#'
#' @slot Rmat A 3-column matrix where each row represents a relative position of
#' interaction. First and second columns represent the position and third column
#' the interaction type.
#' @slot n_neis The number of intracting positions.
#' @slot n_types The number of types of interaction.
setClass("mrfi",
         representation(Rmat = "matrix",
                        n_neis = "numeric",
                        n_types = "numeric"))

setMethod("show", "mrfi",
          function(object){
            cat(object@n_neis, "interacting positions of", object@n_types, "types: \n")
            cat(" rx  | ry  | type \n")
            for(i in seq_len(min(5,object@n_neis))){
              cat(" ",object@Rmat[i,1],"   ", object@Rmat[i,2], "    ", object@Rmat[i,3], "\n")
            }
            if(object@n_neis > 5) {
              cat("... and", (object@n_neis-5), "more.")
            }
          })
