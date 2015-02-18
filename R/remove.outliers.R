

remove.outliers<-function(x, data, cut=3.29){
  
  
  
        outlierMax <- mean(data[,x], na.rm=TRUE) + cut*sd(data[,x], na.rm=TRUE)
        outlierMin <- mean(data[,x], na.rm=TRUE) - cut*sd(data[,x], na.rm=TRUE)
          
        data[,x][data[,x] > outlierMax] <- NA
        data[,x][data[,x] < outlierMin] <- NA
          
          
          
  
  
  #outlierMax <- mean(data[,x], na.rm=TRUE) + 3.29*sd(data[,x], na.rm=TRUE)
  #outlierMin <- mean(data[,x], na.rm=TRUE) - 3.29*sd(data[,x], na.rm=TRUE)
  #n <- length(data[,x])
  #for (i in 1:n){
  #  if(is.na(data[i,x]) == TRUE){
  #    data[i,x] = NA
      
  #  }
    
  #  else if(data[i, x] > outlierMax){
  #    data[i, x] = NA
  #  }
  #  else if (data[i, x] < outlierMin){
  #    data[i, x] = NA
  #  }
  #}  
  return(data)
}

#data = read.csv(file=file.choose())
#data <- remove.outliers("Outliers1", data)




