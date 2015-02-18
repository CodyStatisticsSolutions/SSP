

scale.to.nominal <- function(variable, value, label, data){
  
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


#data <- read.csv(file=file.choose())


#data2 <- scale.to.nominal("Paired1", 20, "Red", data)
#data2 <- scale.to.nominal("Paired1", 19, "Blue", data2)
#dataset.info(data2)

#data2 <- scale.to.nominal("Paired.1", 20, "Red", data)
