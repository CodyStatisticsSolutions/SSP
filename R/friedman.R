friedman <- function(dvs, id, data){
  
  #library(mvtnorm)
  
  dataframe <- cbind(data[id], data[dvs])
    
  myexpr.df <- make.rm(constant=c(id), repeated=c(dvs), data=dataframe);
  
  ranks <- rank(myexpr.df$repdat);
  myexpr.df <- cbind(myexpr.df, ranks);
  means <- lapply(split(myexpr.df$ranks, myexpr.df$contrasts), mean, na.rm=TRUE);
  
  friedman <- lme(ranks~contrasts, data=myexpr.df, random=~1|cons.col);
  sfriedman <- summary(friedman)
  anova <- anova(friedman)
  
  pairwise <- summary(glht(friedman, linfct = mcp(contrasts = "Tukey")))
  
  output <- list()
  
  output$anova2 <- as.list(unclass(anova));
  output$pairwise <- as.list(unclass(pairwise$test$coefficients));
  output$pairwise2 <- as.list(unclass(pairwise$test$pvalues));
  output$means <- as.list(means)
  return(output)
  
}

#library(nlme)
#library(multcomp)
#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#library(opencpu.encode)


#friedman(dvs = c("Paired1", "Paired2Different", "Paired2NotDifferent"), id="ID", data=data);
#friedman(dvs= c("Paired1", "Paired2NotDifferent"), id="ID", data=data)