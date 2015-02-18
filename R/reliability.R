
reliability <- function(vars, data){

	scale <- cbind(data[vars])

	alpha <- alpha(scale);
	items <- length(alpha$keys)
	alpha <- alpha$total[1]

	
	output <- list()
	output$alpha <- unclass(alpha)
	output$items <- items

	return(output);
}
#library(psych)

#reliability(c("Paired1", "Paired2Different", "Paired2NotDifferent"), data)

#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")

