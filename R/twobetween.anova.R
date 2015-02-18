

twobetween.anova <- function(x1, x2, y, data, ...){
  #remove NA
  data <- na.omit(data[c(x1, x2 ,y)]);
  
  #calculations
  myformula <- as.formula(print(paste(y, "~", x1, "+", x2, "+", x1,":", x2)));
  formula2 <- as.formula(paste(y, "~", x1, "*", x2))
  mylm <- lm(myformula, data=data, ...);
  myLevene <- leveneTest(formula2, data=data, center=median);
  myShapiro <- shapiro.test(data[[y]]);

  standard.deviations1 <- lapply(split(data[[y]], data[[x1]]), sd,na.rm=TRUE);
  standard.deviations2 <- lapply(split(data[[y]], data[[x2]]), sd,na.rm=TRUE);
  standard.deviations3 <- lapply(split(data[[y]], data[[x1]]:data[[x2]]), sd,na.rm=TRUE);  
  
  
  means1 <- lapply(split(data[[y]], data[[x1]]), mean, na.rm=TRUE);
  means2 <- lapply(split(data[[y]], data[[x2]]), mean, na.rm=TRUE);
  means3 <- lapply(split(data[[y]], data[[x1]]:data[[x2]]), mean, na.rm=TRUE);  
  
  
  n.length1 <- lapply(split(data[[y]], data[[x1]]), length);
  n.length2 <- lapply(split(data[[y]], data[[x2]]), length);
  n.length3 <- lapply(split(data[[y]], data[[x1]]:data[[x2]]), length);  
  
  
  #boxplot(myformula, data=data);

  myanova <- as.data.frame(Anova(mylm, type=3, singular.ok=T));
  output <- list()

  myanova.ss <- myanova$"Sum Sq";
  Eta <- myanova.ss/(myanova.ss + myanova.ss[length(myanova.ss)])
  myanova$"Mean Sq" <- myanova$"Sum Sq" / myanova$"Df"
  myanova <- myanova[,c(2, 1, 5, 3, 4)]
  myanova <- myanova[-c(1),]

  myanova.ss <- myanova$"Sum Sq";
  Eta <- myanova.ss/(myanova.ss + myanova.ss[length(myanova.ss)])
  myanova$"Mean Sq" <- myanova$"Sum Sq" / myanova$"Df"
  #anova output
  output <- list();
  colnames(myanova) <- c("DF", "SumSq", "MeanSq", "F-value", "P");
  means1 <- lapply(means1, unname);
  means2 <- lapply(means2, unname);
  means3 <- lapply(means3, unname);
  standard.deviations1 <- lapply(standard.deviations1, unname);
  standard.deviations2 <- lapply(standard.deviations2, unname);
  standard.deviations3 <- lapply(standard.deviations3, unname);
  n.length1 <- lapply(n.length1, unname);
  n.length2 <- lapply(n.length2, unname);
  n.length3 <- lapply(n.length3, unname);
  pairwise1 <- pairwise.t.test(data[[y]], data[[x1]], p.adj="bonf");
  pairwise2 <- pairwise.t.test(data[[y]], data[[x2]], p.adj="bonf");
  pairwise3 <- pairwise.t.test(data[[y]], data[[x1]]:data[[x2]], p.adj="bonf");
  output$terms <- apply(myanova,1,as.list);
  output$Eta <- c(Eta[1], Eta[2], Eta[3])  
  output$levene <- myLevene$`Pr(>F)`[1];
  output$shapiro <- myShapiro$p.value;
  output$pairwise1 <- pairwise1$p.value;
  output$pairwise2 <- pairwise2$p.value;
  output$pairwise3 <- pairwise3$p.value;
  output$means1 <- means1;
  output$means2 <- means2;
  output$means3 <- means3;
  output$standard.deviations1 <- standard.deviations1;
  output$standard.deviations2 <- standard.deviations2;
  output$standard.deviations3 <- standard.deviations3;
  output$n1 <- n.length1;
  output$n2 <- n.length2;
  output$n3 <- n.length3;

  
  return(output);
}

#library(car)
#library(foreign)

#data <- read.csv(file=file.choose())


#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#twobetween.anova("Diagnostician", "Dx", "Clean", data)

#twobetween.anova("Not.Different.Group.2", "Chi.Group.1", "No.Correlate", data)
#twobetween.anova("Gender", "Condition", "Clean", data2)

