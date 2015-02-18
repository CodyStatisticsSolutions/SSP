

mann.whitney <- function (x, y, data, ...){

	
	formula <- as.formula(paste(y, "~", x));
	test <- wilcox.test(formula, paired=FALSE, data=data, correct=FALSE);

		
	r <- rank(data[[y]]);
	data2 <- data.frame(r, data[[x]]);
	ranks <- lapply(split(r, data[[x]]), mean, na.rm=TRUE);
	summary <- lapply(split(data[[y]], data[[x]]), summary)

	

	boxplot(formula, data = data, main=paste("Scores on", y, "by", x), xlab = x, ylab = y);

	z <-qnorm(test$p.value/2);

	output <- unclass(test);
	output$ranks <- lapply(ranks, unname);
	output$summary <- lapply(summary, unname);
	output$z <- z;
	return(output);
	

}

#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#mann.whitney("DifferentGroup2", "Normal", data)
