#' Levene test
#' 
#' Perform a levene test
#' 
#' @param y a formula
#' @param data a data frame
#' @return levene test object 
#' @examples levene.test(Sepal.Width ~ Species, data=iris)
#' @import car
#' @export
levene.test <- function(y, data){
	myLevene <- leveneTest(y=y, data=data, center=mean);
	
	output <- list();
	output$formula <- deparse(y);
	output$results <- unclass(myLevene);
	
	return(output);
}
