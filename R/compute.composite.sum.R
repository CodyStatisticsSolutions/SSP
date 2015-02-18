#library(dplyr)

compute.composite.sum <- function(vars, name, data){
  
  variables <- cbind(data[vars])
  data <- mutate(data, newvar = rowSums(variables, na.rm=TRUE))
  vars <- length(colnames(data))
  colnames(data)[[vars]] <- paste(name)

  return(data)
}


#data = read.csv(file=file.choose())
#data<-compute.composite.sum(c("Paired2", "Paired3", "Paired4"), "SuperPaired", data)
