#' Pretty breaks function
#'
#' Makes pretty y limits for regular axes
#'
#' @param max Max y value found in data
#' @return Vector c(ymax, number of breaks)
#' @export
prettybreaks <- function(max){

  #checking for 3 breaks

  #trying to make axis just by rounding
  seg3 <- signif(max/3, 1)
  #print(seg3)
  top3 <- 3*seg3
  #print(top3)
  #checking whether the new limit is high enough and otherwise adding more
  while(0.8*top3 < max){
    #print("add more")
    top3 <- top3 + 3*(10^(floor(log10(max/3))))
    #print(top3)
  }

  #checking for 4 breaks

  #trying to make axis just by rounding
  seg4 <- signif(max/4, 1)
  #print(seg4)
  top4 <- 4*seg4
  #print(top4)
  #checking whether the new limit is high enough and otherwise adding more
  while(0.8*top4 < max){
    #print("add more")
    top4 <- top4 + 4*(10^(floor(log10(max/4))))
    #print(top4)
  }


  #checking for 5 breaks

  #trying to make axis just by rounding
  seg5 <- signif(max/5, 1)
  #print(seg5)
  top5 <- 5*seg5
  #print(top5)
  #checking whether the new limit is high enough and otherwise adding more
  while(0.8*top5 < max){
    #print("add more")
    top5 <- top5 + 5*(10^(floor(log10(max/5))))
    #print(top5)
  }

  if(top4 < top5 & top4 <= top3){
    return(c(top4, 4))
  }
  else if(top5 <= top3){
    return(c(top5, 5))
  }
  else(return(c(top3, 3)))

}
