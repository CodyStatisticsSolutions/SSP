

one.sample.t <- function(y, mu, data){
  normality <- shapiro.test(data[[y]])
  normality.p <- unclass(normality)$p.value;
  
  mean <- mean(data[[y]], na.rm=TRUE)
  sd <- sd(data[[y]], na.rm=TRUE)
  
  test <- t.test(data[[y]], mu=mu, data=data)
  t.stat <- as.list(test$statistic)
  df <- as.list(test$parameter)
  p.value <- test$p.value
  
  output <- list()
  output$normality <- normality.p
  output$mean <- mean
  output$sd <- sd
  output$t.stat <- t.stat$t
  output$df <- df$df
  output$p.value <- p.value
  return(output)  
  
}


#one.sample.t("Correlate", 45, data)

#one.sample.t("Normal", 65, data)







#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
