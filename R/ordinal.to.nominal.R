

ordinal.to.nominal <- function(variable, data){
  
  data[,variable] <- factor(data[,variable], ordered=FALSE)  
  
  return(data)
  
}


