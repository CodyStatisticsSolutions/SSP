stepwise.regression <- function(formula, data, method,...){

  data.r <- na.omit(data[all.vars(formula)])
  
  #fit the model
  
  mylm <- do.call("lm",list(formula, data.r, x=TRUE, y=TRUE, ...));
  
  mylm.r <- step(mylm, direction=method)

  
  #get standardized oefficients
  sd.x <- apply(mylm.r$x, 2, sd);
  sd.y <- sd(mylm.r$y);
  std.coef <- coef(mylm.r) * (sd.x / sd.y);	
  
  #some more plots
  print(plot(mylm.r));
  
  #return R object with relevant data
  output <- list();
  output$ANOVA.table <- apply(as.data.frame(anova(mylm.r)),1,as.list);
  output$COEF.table <- apply(cbind(std.coef,as.data.frame(summary(mylm.r)$coefficients)), 1, as.list);

  
  #model test statistics
  Fstat <- as.list(summary(mylm.r)$fstatistic);
  names(Fstat) <- c("F", "df.numerator", "df.denominator")
  output$model <- Fstat;
  output$model$r.squared <- summary(mylm.r)$r.squared;
  output$model$adj.r.squared <- summary(mylm.r)$adj.r.squared;
  output$model$residual.standard.error <- summary(mylm.r)$sigma;
  output$model$p.value <- pf(Fstat[[1]],Fstat[[2]],Fstat[[3]], lower.tail=FALSE)

  
  
  return(output);
}



#data<- read.csv(file=file.choose())
#library(MASS)
#stepwise.regression(Paired3~Group3 + Paired1 + Paired2Different + Correlate, data, method="backward")
#stepwise.regression(Paired3~Group3 + Paired1 + Paired2Different + Correlate, data, method="both")
#stepwise.regression(Paired3~Group3 + Paired1 + Paired2Different + Correlate, data, method="forward")
#stepwise.regression(NormalDependent~Paired1 + Paired2Different + Correlate, data, method="backward")

#all.vars(NormalDependent~Group3 + Paired1 + Paired2Different + Correlate)
#data3 <- na.omit(data[all.vars(NormalDependent~Group3 + Paired1 + Paired2Different + Correlate)])

#model1 <- lm(NormalDependent~Group3 + Paired1 + Paired2Different + Correlate, data3)
#model1.step <- step(model1, method="backwards")

#model1.step

#stepwise.regression(NormalDependent~Group3 + Paired1 + Paired2Different + Correlate, data, "both")
