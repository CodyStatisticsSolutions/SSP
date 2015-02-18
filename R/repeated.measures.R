repeated.measures <- function(dvs, id, data){
  
  #library(mvtnorm)
	
	dataframe <- cbind(data[id], data[dvs])

	mean <- sapply(dataframe, mean, na.rm=TRUE)
	sd <- sapply(dataframe, sd, na.rm=TRUE)

	shapiros <- sapply(data[dvs], shapiro.test);

	myexpr.df <- make.rm(constant=c(id), repeated=c(dvs), data=dataframe);
	rmanova <- lme(repdat~contrasts, data=myexpr.df, random=~1|cons.col, na.action=na.omit);
	srmanova <- summary(rmanova)
	anova <- anova(rmanova)

	pairwise <- summary(glht(rmanova, linfct = mcp(contrasts = "Tukey")))

	output <- list()

	output$mean <- as.list(unclass(mean));
	output$sd <- as.list(unclass(sd));
	output$shapiros <- as.list(as.data.frame(shapiros));
	output$anova2 <- anova;
	output$pairwise <- as.list(unclass(pairwise$test$coefficients));
	output$pairwise2 <- as.list(unclass(pairwise$test$pvalues));
	return(output)
	
}

#library(nlme)
#library(multcomp)

#data2 <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#library(opencpu.encode)


#repeated.measures(dvs = c("Paired1", "Paired2Different", "Paired2NotDifferent"), id="ID", data=data2);

#cat(asJSON(repeated.measures(dvs = c("Paired1", "Paired2Different", "Paired2NotDifferent"), id="ID", data=data)));