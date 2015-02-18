
change.to.numeric <- function(variable, data){
  
  data[,variable] <- as.numeric(data[,variable])  
  
  return(data)
  
}