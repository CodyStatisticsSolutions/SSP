

scale.to.ordinal <- function(variable, value, label, data){
  
  n <- length(data[,variable])
  
  data[,variable] <- as.character(data[,variable])
  
  for (i in 1:n){
    if(is.na(data[i, variable])){
      data[i, variable] <- NA
    }
    else if(data[i, variable] == value){
      data[i, variable] <- label
    }
  }  
  
  data[,variable] <- as.factor(data[,variable])
  return(data)
}

#data2 <- scale.to.ordinal("TypicalLikert1", 1, "A", data)
#data2 <- scale.to.ordinal("TypicalLikert1", 2, "B", data2)
#data2 <- scale.to.ordinal("TypicalLikert1", 3, "C", data2)
#data2 <- scale.to.ordinal("TypicalLikert1", 4, "D", data2)
#data2 <- scale.to.ordinal("TypicalLikert1", 5, "E", data2)

#data2 <- nominal.to.ordinal("TypicalLikert1", c("A", "B", "C", "D", "E"), data2)
#data2$TypicalLikert1
