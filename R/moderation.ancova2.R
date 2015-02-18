
moderation.ancova2 <- function(iv, mv, dv, data){
	
	formula <- as.formula(paste(dv, "~as.factor(", iv, ") +",mv,"+ as.factor(",iv,")*",mv));

	mylm <- lm(formula, data=data,x=TRUE, y=TRUE);

	
	#get standardized oefficients
	#sd.x <- apply(mylm$x, 2, sd);
	#sd.y <- sd(mylm$y);
	#std.coef <- coef(mylm) * (sd.x / sd.y);


	#return R object with relevant data
	output <- list();
  
  myanova <- as.data.frame(Anova(mylm, type = 3))
	myanova$"Mean Sq" <- myanova$"Sum Sq" / myanova$"Df"
	myanova <- myanova[,c(2, 1, 5, 3, 4)]
	myanova <- myanova[-c(1),]  
  

  output$ANOVA.table <- apply(myanova,1, as.list)
	
  
  
  #output$ANOVA.table2 <- apply(as.data.frame(anova(mylm)),1,as.list);
	#output$COEF.table <- apply(cbind(std.coef,as.data.frame(summary(mylm)$coefficients)), 1, as.list);


	#model test statistics
	Fstat <- as.list(summary(mylm)$fstatistic);
	names(Fstat) <- c("F", "df.numerator", "df.denominator")
	output$model <- Fstat;
	#output$model$r.squared <- summary(mylm)$r.squared;
	#output$model$adj.r.squared <- summary(mylm)$adj.r.squared;
	#output$model$residual.standard.error <- summary(mylm)$sigma;
	output$model$p.value <- pf(Fstat[[1]],Fstat[[2]],Fstat[[3]], lower.tail=FALSE)

	
	return(output);
}

#moderation.ancova2("Group.3", "IVModerate", "Not.Normal.Dependent", data)
#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#data2 <- read.csv(file=file.choose())

