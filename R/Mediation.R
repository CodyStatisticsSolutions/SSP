


mediation <- function(iv, mv, dv, data){

	formula1 <- as.formula(paste(dv,"~",iv));
	formula2 <- as.formula(paste(mv,"~",iv));
	formula3 <- as.formula(paste(dv,"~",iv,"+",mv));
	
	linear1 <- lm(formula1, data=data, x=TRUE, y=TRUE);
	linear2 <- lm(formula2, data=data, x=TRUE, y=TRUE);
	linear3 <- lm(formula3, data=data, x=TRUE, y=TRUE);
	Fstat1 <- as.list(summary(linear1)$fstatistic);
	Fstat2 <- as.list(summary(linear2)$fstatistic);
	Fstat3 <- as.list(summary(linear3)$fstatistic);

	sd.x1 <- apply(linear1$x, 2, sd);
	sd.y1 <- sd(linear1$y);
	std.coef1 <- coef(linear1) * (sd.x1 / sd.y1);


	output<-list();
	output$ANOVA.table1 <- apply(as.data.frame(anova(linear1)),1,as.list);
	output$COEF.table1 <- apply(cbind(std.coef1,as.data.frame(summary(linear1)$coefficients)), 1, as.list);

	sd.x2 <- apply(linear2$x, 2, sd);
	sd.y2 <- sd(linear2$y);
	std.coef2 <- coef(linear2) * (sd.x2 / sd.y2);

	output$ANOVA.table2 <- apply(as.data.frame(anova(linear2)),1,as.list);
	output$COEF.table2 <- apply(cbind(std.coef2,as.data.frame(summary(linear2)$coefficients)), 1, as.list);

	sd.x3 <- apply(linear3$x, 2, sd);
	sd.y3 <- sd(linear3$y);
	std.coef3 <- coef(linear3) * (sd.x3 / sd.y3);

	output$ANOVA.table3 <- apply(as.data.frame(anova(linear3)),1,as.list);
	output$COEF.table3 <- apply(cbind(std.coef3,as.data.frame(summary(linear3)$coefficients)), 1, as.list);



	#model test statistics
	names(Fstat1) <- c("F1", "df.numerator1", "df.denominator1")
	output$model1 <- Fstat1;
	output$model1$r.squared <- summary(linear1)$r.squared;
	output$model1$adj.r.squared <- summary(linear1)$adj.r.squared;
	output$model1$residual.standard.error <- summary(linear1)$sigma;
	output$model1$p.value <- pf(Fstat1[[1]],Fstat1[[2]],Fstat1[[3]], lower.tail=FALSE)

	names(Fstat2) <- c("F2", "df.numerator2", "df.denominator2")
	output$model2 <- Fstat2;
	output$model2$r.squared <- summary(linear2)$r.squared;
	output$model2$adj.r.squared <- summary(linear2)$adj.r.squared;
	output$model2$residual.standard.error <- summary(linear2)$sigma;
	output$model2$p.value <- pf(Fstat2[[1]],Fstat2[[2]],Fstat2[[3]], lower.tail=FALSE)

	names(Fstat3) <- c("F3", "df.numerator3", "df.denominator3")
	output$model3 <- Fstat3;
	output$model3$r.squared <- summary(linear3)$r.squared;
	output$model3$adj.r.squared <- summary(linear3)$adj.r.squared;
	output$model3$residual.standard.error <- summary(linear3)$sigma;
	output$model3$p.value <- pf(Fstat2[[1]],Fstat2[[2]],Fstat2[[3]], lower.tail=FALSE)
	output$model3$p.value <- pf(Fstat3[[1]],Fstat3[[2]],Fstat3[[3]], lower.tail=FALSE)
	

	return(output);



	
}


#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#mediation("Correlate", "NormalDependent","Paired1", data=data)



#library(opencpu.encode);
#cat(asJSON(mediation(Correlate, NormalDependent,Paired1, data=data)));


	