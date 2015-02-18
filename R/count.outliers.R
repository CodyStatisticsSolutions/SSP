#data = read.csv(file=file.choose())
#data2 <- remove.outliers("Outliers1", data)


count.outliers <- function(x, data){
  outlierMax <- mean(data[,x], na.rm=TRUE) + 3.29*sd(data[,x], na.rm=TRUE)
  outlierMin <- mean(data[,x], na.rm=TRUE) - 3.29*sd(data[,x], na.rm=TRUE)
  n <- length(data[,x])
  
  outlier <- c()
  
  for (i in 1:n){
    if(is.na(data[i,x]) == TRUE){
      outlier[i] = 0
      
    }
    else if(data[i, x] > outlierMax){
      outlier[i] = 1
    }
    else if (data[i, x] < outlierMin){
      outlier[i] = 1
    }
    else{
      outlier[i] = 0
    }
  }
  return(sum(outlier, na.rm=TRUE))
  
}

#data = read.csv(file=file.choose())
#count.outliers("Paired1", data)
