#' Student t-test
#' 
#' Compare two groups in a t-test.
#' 
#' @param x factor variable with exactly two levels
#' @param y numeric variable
#' @param data data frame
#' @param ... arguments pased on to t.test
#' @return bunch of t-test info
#' @import stats
#' @examples student.test("vs", "mpg", data=mtcars)
#' @export
student.test <- function(x, y, data, ...){
	## Note that x is a factor variable with exactly 2 levels!
	myformula <- as.formula(paste(y, "~", x));
	output <- t.test(myformula, data, ...);
	boxplot(myformula, data);

	standard.deviations <- lapply(split(data[[y]], data[[x]]), sd,na.rm=TRUE);
	means <- lapply(split(data[[y]], data[[x]]), mean,na.rm=TRUE);
	sd.pooled <- sqrt(mean(unlist(standard.deviations)^2));
	cohen.d <- abs(diff(unlist(means))) / sd.pooled;

	output <- unclass(output);	
	standard.deviations <- lapply(standard.deviations, unname);
	means <- lapply(means, unname);
	output$standard.deviations <- standard.deviations;
	output$means <- means;
	output$cohen.d <- unname(cohen.d);
	
	#shapiro test over full y variable
	output$shapiro.p.value <- shapiro.test(data[[y]])$p.value
	output$shapiro.groups <- shapiro_internal(x, y, data)
	output$t.test2 <- t.test(myformula, data, var.equal=TRUE, ...);
	output$levenes.p.value <- leveneTest(data[[y]] ~ as.factor(data[[x]]), center=mean)$`Pr(>F)`
	output$t.test2 <- unclass(output$t.test2);
	return(output);
}


#library(StatSolutions)
#library(car)
#student.test("vs", "mpg", data=mtcars)
#lev <- leveneTest(mtcars$mpg ~ as.factor(mtcars$vs), data=mtcars)
#unclass(lev)
#library(opencpu.encode);
#cat(asJSON(student.test("vs", "mpg", data=mtcars)));
