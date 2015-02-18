#' Linear Regression
#' 
#' Perform a linear regression analysis with one or more predictors
#' 
#' @param formula regression formula
#' @param data data frame
#' @param ... other stuff passed to lm()
#' @return lots of output
#' @examples linear.regression(speed~dist, data=cars)
#' @export
linear.regression <- function(formula, data, ...){
	
	#fit the model
	mylm <- lm(formula, data, x=TRUE, y=TRUE, ...);
  #vifs <- vif(mylm)
	
	#get standardized oefficients
	sd.x <- apply(mylm$x, 2, sd);
	sd.y <- sd(mylm$y);
	std.coef <- coef(mylm) * (sd.x / sd.y);	
	
	#some more plots
	print(plot(mylm));
	
	#return R object with relevant data
	output <- list();
	output$ANOVA.table <- apply(as.data.frame(anova(mylm)),1,as.list);
	output$COEF.table <- apply(cbind(std.coef,as.data.frame(summary(mylm)$coefficients)), 1, as.list);
	
	#model test statistics
	Fstat <- as.list(summary(mylm)$fstatistic);
	names(Fstat) <- c("F", "df.numerator", "df.denominator")
	output$model <- Fstat;
	output$model$r.squared <- summary(mylm)$r.squared;
	output$model$adj.r.squared <- summary(mylm)$adj.r.squared;
	output$model$residual.standard.error <- summary(mylm)$sigma;
	output$model$p.value <- pf(Fstat[[1]],Fstat[[2]],Fstat[[3]], lower.tail=FALSE)
  #output$vifs <- as.list(vifs)

	
	return(output);
}
#data <-read.csv(file=file.choose())

#linear.regression(NormalDependent~Group3 + Paired1 + Paired2Different + Correlate, data)
#library(car)

#linear.regression(Normal ~ Paired1 + Paired2Different + Correlate, data=data)

