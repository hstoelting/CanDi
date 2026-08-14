#' Pretty breaks (log) function
#'
#' Makes pretty y limits for log axes
#'
#' @param min Min y value found in data (excl. 0)
#' @param max Max y value found in data
#' @return Numeric limit vector c(lower, upper)
#' @export
prettybreaks_log <- function(min, max){

  lower <- 10^(log(min, 10) %>% floor())
  upper <- 10*(10^(log(max, 10) %>% ceiling()))
  # print(lower)
  # print(upper)

  while(0.8*(log(upper, 10)-log(lower, 10)) < (log(max,10)-log(lower, 10))){
    upper <- upper*10
  }

    while(upper/lower < 1000){
      lower <- lower/10
      upper <- upper*10

      # print(lower)
      # print(upper)
    }

  return(c(lower, upper))
}

