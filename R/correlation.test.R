#' Correlation test
#' 
#' Perform a correlation test between two numeric variables
#' 
#' @param x a numeric variable 
#' @param y anoter numeric variable
#' @param data a dataframe
#' @param ... passed on to stats::cor.test
#' @return cor.test object
#' @import stats
#' @export
#' @examples correlation.test("vs", "mpg", data=mtcars)
#' 
correlation.test <- function(x, y, data, ...){
	myformula <- as.formula(paste("~", x, "+", y));
	plot(data[c(x,y)])
	output <- stats::cor.test(myformula, data=data, ...);
	return(unclass(output));	
}
