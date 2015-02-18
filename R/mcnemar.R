mcnemar <- function (x1, x2, data){
  
  mcnemar <- mcnemar.test(data[[x1]],data[[x2]])
  table <- table(data[[x1]], data[[x2]])
  
  table <-as.data.frame(table)
  
  factors1 <- levels(data[[x1]])
  factors2 <- levels(data[[x2]])
  
  repeated <- if(factors1[1] == factors2[1] & factors1[2] == factors2[2]) "met" else "not met"
  
  chi <- mcnemar$statistic
  df <- mcnemar$parameter
  p <- mcnemar$p.value
  
  output <- list()
  output$chi.square <- unname(chi)
  output$p <- p
  output$table <- table$Freq
  output$factor <- factors1
  output$repeated <- repeated
  
  return(output)
}

#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")

#mcnemar("DifferentGroup2", "ChiGroup1", data)

#mcnemar("DifferentGroup2", "NotDifferentGroup2", data)

#mcnemar("DifferentGroup2", "Group2", data)