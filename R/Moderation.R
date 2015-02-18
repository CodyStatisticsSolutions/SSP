
moderation <- function(iv, mv, dv, data){

	ivd <- data[[iv]];
	ivc <- (ivd - mean(ivd));
	mvd <- data[[mv]];
	mvc <- (mvd - mean(mvd));
	interaction <- (ivc*mvc);	
	
	formula <- as.formula(paste(dv,"~",iv,"+",mv,"+","interaction"));

	mylm <- lm(formula, data=data,x=TRUE, y=TRUE, na.rm=T);

	#get standardized coefficients
	sd.x <- apply(mylm$x, 2, sd);
	sd.y <- sd(mylm$y);
	std.coef <- coef(mylm) * (sd.x / sd.y);


	#return R object with relevant data
	output <- list();
	output$ANOVA.table <- apply(as.data.frame(anova(mylm, type="III")),1,as.list);
	output$COEF.table <- apply(cbind(std.coef,as.data.frame(summary(mylm)$coefficients)), 1, as.list);


	#model test statistics
	Fstat <- as.list(summary(mylm)$fstatistic);
	names(Fstat) <- c("F", "df.numerator", "df.denominator")
	output$model <- Fstat;
	output$model$r.squared <- summary(mylm)$r.squared;
	output$model$adj.r.squared <- summary(mylm)$adj.r.squared;
	output$model$residual.standard.error <- summary(mylm)$sigma;
	output$model$p.value <- pf(Fstat[[1]],Fstat[[2]],Fstat[[3]], lower.tail=FALSE)

	
	return(output);
}

#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#attach(data)
#moderation("Correlate", "NormalDependent","Paired1", data=data)



#library(opencpu.encode);
#cat(asJSON(moderation(Correlate, NormalDependent,Paired1, data=data)));

