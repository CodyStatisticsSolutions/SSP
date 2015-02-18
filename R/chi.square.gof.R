chi.square.gof <- function(x, data){
  xvar <- as.factor(data[[x]]);
  
  chitest <- unclass(chisq.test(table(xvar)));
  chitest <- lapply(chitest, unclass); #unclass all individual elements as well

  
  
  chitest$dimnames <- dimnames(chitest$observed);
  chitest$observed <- unclass(chitest$observed);
  chitest$expected <- unclass(chitest$expected);
  
  
  return(chitest);
}



#chi.square.gof("ChiGroup2Large", data=data)


