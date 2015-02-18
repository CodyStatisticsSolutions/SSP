

#library(psy)


cronbach.alpha <- function(vars, newvar, data){
  variables <- cbind(data[vars])
  cronbach <- cronbach(variables)
  
  mean <- mean(data[,newvar], na.rm=TRUE)
  sd <- sd(data[,newvar], na.rm=TRUE)
  
  output <- list()
  output$cronbach <- cronbach
  output$mean <- mean
  output$sd <- sd
  return(output)  
}

#data = read.csv(file=file.choose())
#cronbach.alpha(c("TypicalLikert1.R", "TypicalLikert2", "TypicalLikert3.R"), "TypicalComposite", data)