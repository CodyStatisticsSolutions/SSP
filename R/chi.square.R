#' Chi Square test
#' 
#' Perform a Chi Square test
#'
#' @param x A string defining a categorial variable that exists in the dataset. 
#' @param y A string defining another categorial variable that exists in the dataset. 
#' @param data The dataset
#' @import stats
#' @return A valid JSON string
#' @export
#' @examples chi.square("vs", "am", data=mtcars)
#'

chi.square <- function(x, y, data){
	xvar <- as.factor(data[[x]]);
	yvar <- as.factor(data[[y]]);
	
	chitest <- unclass(chisq.test(xvar,yvar,correct=FALSE));
	chitest <- lapply(chitest, unclass); #unclass all individual elements as well
  
  fisher <- fisher.test(xvar, yvar)$p.value
  
	
	chitest$dimnames <- dimnames(chitest$observed);
	chitest$observed <- unclass(chitest$observed);
	chitest$expected <- unclass(chitest$expected);
  chitest$fisher <- fisher

	
	return(chitest);
}

#chi.square("vs", "am", data=mtcars)

#library(opencpu.encode)
#cat(asJSON(chi.square("vs", "am", data=mtcars)))

#chi.square("DifferentGroup2", "ChiGroup2Small", data=data)
#chi.square("Condition", "Gender", data=data3)

