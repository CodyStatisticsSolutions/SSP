#library(dplyr)

compute.composite.mean <- function(vars, name, data){
  
  variables <- cbind(data[vars])
  data <- mutate(data, newvar = rowMeans(variables, na.rm=TRUE))
  vars <- length(colnames(data))
  colnames(data)[[vars]] <- paste(name)
  
  return(data)
}


#data = read.csv(file=file.choose())
#data<-compute.composite.mean(c("TypicalLikert1.R", "TypicalLikert2", "TypicalLikert3.R"), "TypicalComposite", data)