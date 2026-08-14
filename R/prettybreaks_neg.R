#' Pretty breaks (neg) function
#'
#' Makes pretty y limits for regular axes that dip below zero
#'
#' @param min Min y value found in data (less than 0)
#' @param max Max y value found in data
#' @return Vector c(ymax, number of breaks)
#' @export
prettybreaks_neg <- function(min, max){

  #checking for 2+x breaks

  #trying to make axis just by rounding
  seg2 <- signif(max/2, 1)

  top2 <- 2*seg2
  bot2 <- -seg2
  n2 <- 2
  r2 <- top2-bot2

  #checking whether the new limit is high enough and otherwise adding more
  while((0.8*(top2-bot2) < (max-bot2)) | top2 < max | bot2 > min){
    #print("help")
    while(top2 < max){
      top2 <- top2 + 2*(10^(floor(log10(max/2))))
      seg2 <- top2/2
      bot2 <- -seg2
    }

    while(bot2 > min){
      bot2 <- bot2 - seg2

    }

    while((0.8*(top2-bot2) < (max-bot2))){
      top2 <- top2 + 2*(10^(floor(log10(max/2))))
      seg2 <- top2/2
      bot2 <- -seg2
    }

    n2 <- (top2-bot2)/seg2
    r2 <- top2-bot2

  }

  #checking for 3+x breaks

  #trying to make axis just by rounding
  seg3 <- signif(max/3, 1)

  top3 <- 3*seg3
  bot3 <- -seg3
  n3 <- 3
  r3 <- top3-bot3

  #checking whether the new limit is high enough and otherwise adding more
  while((0.8*(top3-bot3) < (max-bot3)) | top3 < max | bot3 > min){
    #print("help")
    while(top3 < max){
      top3 <- top3 + 3*(10^(floor(log10(max/3))))
      seg3 <- top3/3
      bot3 <- -seg3
    }

    while(bot3 > min){
      bot3 <- bot3 - seg3

    }

    while((0.8*(top3-bot3) < (max-bot3))){
      top3 <- top3 + 3*(10^(floor(log10(max/3))))
      seg3 <- top3/3
      bot3 <- -seg3
    }

    n3 <- (top3-bot3)/seg3
    r3 <- top3-bot3

  }

  #checking for 4+x breaks

  #trying to make axis just by rounding
  seg4 <- signif(max/4, 1)

  top4 <- 4*seg4
  bot4 <- -seg4
  n4 <- 4
  r4 <- top4-bot4

  #checking whether the new limit is high enough and otherwise adding more
  while((0.8*(top4-bot4) < (max-bot4)) | top4 < max | bot4 > min){
    #print("help")
    while(top4 < max){
      top4 <- top4 + 4*(10^(floor(log10(max/4))))
      seg4 <- top4/4
      bot4 <- -seg4
    }

    while(bot4 > min){
      bot4 <- bot4 - seg4

    }

    while((0.8*(top4-bot4) < (max-bot4))){
      top4 <- top4 + 4*(10^(floor(log10(max/4))))
      seg4 <- top4/4
      bot4 <- -seg4
    }

    n4 <- (top4-bot4)/seg4
    r4 <- top4-bot4

  }

  if(r3 <= r4 & r3 <= r2){
    return(c(bot3, top3, n3))
  }
  else if(r4 <= r2){
    return(c(bot4, top4, n4))
  }
  else{return(c(bot2, top2, n2))}

}

prettybreaks_neg(min = -3.2, max = 4)
