stack.data <- function(dvs, ivs, id, data){
  
  dataframe <- cbind(data[id], data[dvs], data[ivs])
  
  stacked <- make.rm(constant=c(id, ivs), repeated=c(dvs), data=dataframe);
  
  stacked$name <- stacked$repdata
  
  return(stacked)
}

#library(nlme)
#data <- read.csv(file=file.choose())
#data2 <- stack.data(dvs=c("Paired1", "Paired2", "Paired3", "Paired4"), ivs=c("Group3", "Correlate"), id = "ID", data=data)
#write.csv(data2, file="C:/Users/Statistics Solutions/Desktop/Examples/Stacked.csv")
