nominal.to.ordinal <- function (x, order, data){
  
  
  data[[x]] <- ordered(data[[x]], levels = order)
  return(data)
  
}


#data2 <- nominal.to.ordinal("ChiGroup2Large", order=c("Alpha", "Beta", "Theta"), data)

