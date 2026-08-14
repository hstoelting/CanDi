library(magrittr)
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

prettybreaks_log2 <- function(min, max){
  
  lower <- min %>% floor()
  upper <- max %>% ceiling()
  # print(lower)
  # print(upper)
  
  while(0.8*(upper-lower) < max-lower){
    upper <- upper+1
  }
  
  while(upper-lower < 3){
    lower <- lower-1
    if(upper-lower < 3){
      upper <- upper+1
    }
    
    
    # print(lower)
    # print(upper)
  }
  
  range <- upper-lower 

  
  while(range %% 3 != 0 & range %% 4 != 0 & range %% 5 != 0 ){
    lower <- lower -1 
    if(range %% 3 != 0 & range %% 4 != 0 & range %% 5 != 0 ){
      upper <- upper+1
    }
    range <- upper-lower 
  }
    
  if(range %% 4 == 0){nbreaks <- 4}
  else if(range %% 3 == 0){nbreaks <- 3}
  else if(range %% 5 == 0){nbreaks <- 5}

  return(c(lower, upper, nbreaks))
  
}


