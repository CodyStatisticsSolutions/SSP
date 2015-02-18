#' Multiple correlations
#' 
#' Do an analysis with multiple correlations. 
#' 
#' @param data a data frame
#' @param vars (optional) a vector with the subset of variable names
#' @param type type of correlation. Defaults to "pearson".
#' @return correlation matrix
#' @examples multi.cor(iris)
#' @import Hmisc
#' @export
multi.cor <- function(data, vars, type="pearson"){
	mydf <- data;
	if(missing(vars) || vars == ""){
		cordata <- mydf[sapply(mydf, is.numeric)];
	} else {
		cordata <- mydf[unlist(vars)];
	}
	print(plot(cordata));
	output <- rcorr(as.matrix(cordata), type=type);
	output <- lapply(output, as.data.frame)
	return(unclass(output));	
}
