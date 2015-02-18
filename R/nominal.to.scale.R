

nominal.to.scale <- function(variable, label, value, data){
  
  n <- length(data[,variable])  
  
  data[,variable] <- as.character(data[,variable])
  
  for (i in 1:n){
    if(is.na(data[i, variable])){
      data[i, variable] <- NA
    }
    else if(data[i, variable] == label){
      data[i, variable] <- value
    }
  }  

  return(data)
  
}


#data <- read.csv(file=file.choose())


#data2 <- nominal.to.scale("Group2", "Blue", 1, data)
#data2 <- nominal.to.scale("Group2", "Green", 2, data2)
#data2 <- change.to.numeric("Group2", data2)
#dataset.info(data2)
