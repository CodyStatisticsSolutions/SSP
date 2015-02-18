grand.center <- function(var,name,data) {
  require(plyr)
  
 
  data2 <- mutate(data, var.center = data[[var]] - mean(data[[var]],na.rm=T))
  
  
  num.vars <- length(colnames(data2))
  
  colnames(data2)[[num.vars]] <- paste(name)
  
  return(data2)
  
}

#data3 <- grand.center("Correlate", "Correlate.Grand", data2)

