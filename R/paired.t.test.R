#' Paired t-test
#' 
#' Do a t-test on paired observations
#' 
#' @param y1 variable name 1
#' @param y2 variable name 2
#' @param data data frame
#' @param ... arguments passed to t.test
#' @return lots of output
#' @import stats
#' @examples paired.t.test("Sepal.Width", "Sepal.Length", iris);
#' @export
paired.t.test <- function(y1, y2, data, ...){
	
	mydata <- data[c(y1,y2)];
	boxplot(mydata);
	output <- t.test(mydata[[y1]], mydata[[y2]], paired=T, ...);
	y3 <- mydata[[y2]] - mydata [[y1]]
	norm1 <- shapiro.test(y3);
	output$shapiro.p.value1 <- unclass(norm1)$p.value;
	
	Datas <- data.frame(mydata[[y1]], mydata[[y2]])
	Datas2 <- Datas[complete.cases(Datas),]
	
	output$Mean1 <- mean(Datas2[,1]);
	output$Mean2 <- mean(Datas2[,2]);
	output$SD1 <- sd(Datas2[,1]);
	output$SD2 <- sd(Datas2[,2]);
	Diff <- c(mydata[[y1]]) - c(mydata[[y2]]);
	M.Diff <- mean(Diff, na.rm=TRUE);
	SD.Diff <- sd(Diff, na.rm=TRUE);
	output$Cohen.d <- abs(M.Diff/SD.Diff);
	return(unclass(output));
}

#library(opencpu.encode)
#cat(asJSON(paired.t.test("Sepal.Width", "Sepal.Length", iris)));
#t.test(iris$Sepal.Width,iris$Sepal.Length, data=iris, paired=T)
