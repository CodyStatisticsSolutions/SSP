


ordinal.to.scale <- function(variable, label, value, data){
  
  n <- length(data[,variable])  
  
  data[,variable] <- as.character(data[,variable])
  
  for (i in 1:n){
    if(data[i, variable] == label){
      data[i, variable] <- value
    }
  }  
  
  return(data)
  
}

#data3 <- ordinal.to.scale("TypicalLikert1", "A", 1, data2)
#data3 <- ordinal.to.scale("TypicalLikert1", "B", 2, data3)
#data3 <- ordinal.to.scale("TypicalLikert1", "C", 3, data3)
#data3 <- ordinal.to.scale("TypicalLikert1", "D", 4, data3)
#data3 <- ordinal.to.scale("TypicalLikert1", "E", 5, data3)

#data3 <- change.to.numeric("TypicalLikert1", data3)

#data3$TypicalLikert1