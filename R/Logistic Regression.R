#library(fmsb)

logistic.regression <- function(formula, dv, data){

	myglm <- glm(formula, data=data, family=binomial);
	sum.coef <-(summary(myglm)$coef)

	

	est <- exp(sum.coef[,1])
	upper.ci <- exp(sum.coef[,1] + 1.96*sum.coef[,2])
	lower.ci <- exp(sum.coef[,1] - 1.96*sum.coef[,2])

	chi.square <- (summary(myglm)$null.deviance) - (summary(myglm)$deviance);
	df <- (summary(myglm)$df.null) - (summary(myglm)$df.residual)
	pvalue <- pchisq(chi.square, df, lower.tail=FALSE)
	levels <- levels(data[[dv]])
		
	output <- list()
	output$terms <- lapply(levels, unname)
	output$coefficients <- apply(as.data.frame(summary(myglm)$coefficients),1,as.list);
	output$ci <- apply(as.data.frame(cbind(est, lower.ci, upper.ci)),1,as.list);
	output$chi.square <- chi.square;
	output$df <- df;
	output$pvalue <- pvalue;
	output$nagelkerke <- NagelkerkeR2(myglm)$R2;
	return(unclass(output));
	
}

#logistic.regression(DifferentGroup2 ~ NotNormalDependent + ID + as.factor(Group3), "DifferentGroup2", data=data)
#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")


#logistic.regression(NotDifferentGroup2 ~ Paired1 + as.factor(ChiGroup1), "NotDifferentGroup2", data=data)