
ancova <- function (x, cv, y, data){
  
  require(mvtnorm)
  require(multcomp)
  require(effects)
  require(car)
	
	mydata <- na.omit(data[c(x,cv,y)]);

	formula <- as.formula(paste(y,"~",x,"+",cv));
	formula.2 <- as.formula(paste(y,"~",x));
	mylm <- lm(formula, data=mydata);
	myLevene <- leveneTest(formula.2, data=mydata, center=median);
	myShapiro <- shapiro.test(mydata[[y]]);
	means <- lapply(split(mydata[[y]], mydata[[x]]), mean, na.rm=TRUE);
	sds <- lapply(split(mydata[[y]], mydata[[x]]), sd, na.rm=TRUE);
	n.length <- lapply(split(mydata[[y]], mydata[[x]]), length);

	myanova <- as.data.frame(anova(mylm));
	
	myanova.ss <- myanova$"Sum Sq";
	Eta <- myanova.ss/(myanova.ss + myanova.ss[length(myanova.ss)])
        
        print(qqnorm(mylm$res))
        print(qqline(mylm$res))
        
        

	output <- list();
	colnames(myanova) <- c("DF", "SumSq", "MeanSq", "F-Value", "p");
	means <- lapply(means, unname);
	sds <- lapply(sds, unname);
	n.length <- lapply(n.length, unname);
  
  adjmeans <- effect(x, mylm, se=TRUE)
  
  
  

	linfcts <- eval(parse(text=paste("mcp(",paste(x,' ="Tukey" ',collapse=','),')')));
	pairwise <- summary(glht(mylm , linfct = linfcts));

	pairwise2 <- as.list(unclass(pairwise$test$coefficients));
	pairwise3 <- as.list(unclass(pairwise$test$pvalues));

	pairwise4 <- rbind(pairwise2, pairwise3)

	output$terms <- apply(myanova,1,as.list);
	output$Eta1 <- Eta[1]
	output$Eta2 <- Eta[2]
	output$levene <- myLevene$`Pr(>F)`[1];
	output$shapiro <- myShapiro$p.value;
	output$pairwise <- as.list(unclass(pairwise4[2,]));
	output$means <- means;
	output$sds <- sds;
	output$n <- n.length;
  output$adjmeans <- as.list(unname(adjmeans$fit))
	return(output);
}


#library(car)
#library(foreign)
#library(multcomp)
#library(effects)


#data <- read.csv(file.choose())


#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")

#ancova("Group3", "Paired1", "NotNormalDependent", data)
#ancova("ChiGroup2Large", "Paired1", "Correlate", data)

#library(jsonlite)
#toJSON(ancova("Group3", "Paired1", "NotNormalDependent", data))


