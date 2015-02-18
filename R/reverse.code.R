
#library(dplyr)


reverse.code <- function(var, min, max, data){

  data <- mutate(data, newvar = max + min - data[[var]])
  vars <- length(names(data))
  newname <- paste(var, ".R", sep="")
  colnames(data)[[vars]] <- newname
  return(data)
}


#data <- read.csv(file=file.choose())
#data <- reverse.code("TypicalLikert3", 1, 5, data)
