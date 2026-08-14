#' Asterisk function
#'
#' Takes a p value and returns the right number of asterisks
#'
#' @param p p value (numeric)
#' @return ns or correct number of * (character)
#' @export
asterisk <- function(p) {
  if (is.numeric(p) == FALSE){print("p value must be numeric")
    return(NA)}
  else if(is.na(p) == TRUE){return(NA)}
  else if(p < 0.0001) {return("****")}
  else if (p < 0.001) {return ("***")}
  else if (p < 0.01) {return ("**")}
  else if (p < 0.05) {return("*")}
  else {return ("ns")}
}
