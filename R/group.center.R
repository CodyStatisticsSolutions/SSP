group.center <- function(var,grp,name,data) {
  require(plyr)
  
  data <- mutate(data, var.center = data[[var]] - ave(data[[var]], data[[grp]], FUN =function(x) mean(x, na.rm=TRUE)))  
  
  num.vars <- length(colnames(data))
  colnames(data)[[num.vars]] <- paste(name)
  
  return(data)
  
}


#data4 <- group.center("Correlate", "contrasts", "Correlate.Mean", data1)
