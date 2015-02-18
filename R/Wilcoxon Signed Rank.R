
wilcox.signed.rank <- function(y1, y2, data, ...){
	
	y1data <- data[[y1]];
	y2data <- data[[y2]];
	datas <- as.data.frame(c(y1data, y2data))
	data2 <- datas[complete.cases(datas),]

	boxplot(y1data, y2data, names=c(y1, y2), data=data2, ylab="Scores", main=paste("Boxplot for",y1,"and",y2));
	test <- wilcox.test(y1data, y2data, paired=TRUE, data=data2);
	
	med1 <- median(y1data);
	med2 <- median(y2data);
	summary1 <- as.list(summary(y1data));
	summary2 <- as.list(summary(y2data));
	output <- unclass(test);
	output$med1 <- med1;
	output$med2 <- med2;
	output$sum1 <- summary1;
	output$sum2 <- summary2;
	output$z <- qnorm(test$p.value)
	return(output)
	
}


#wilcox.signed.rank("Paired1", "Paired2Different", data=data)
#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")


