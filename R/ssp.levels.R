ssp.levels <- function (variable, data){
  
  level <- levels(as.factor(data[[variable]]))
  return(as.list(level))
}