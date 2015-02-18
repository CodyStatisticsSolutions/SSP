#' Oneway Anova
#' 
#' Analysis of Variance for a linear model
#' 
#' @param x Variable name 1
#' @param y Variable name 2
#' @param data data frame
#' @param ... arguments passed on to lm()
#' @return lots of output
#' @import car
#' @examples oneway.anova("Species", "Sepal.Width", iris)
#' @export
oneway.anova <- function(x, y, data, ...){
	#remove NA
	data <- na.omit(data[c(x,y)]);
	
	#calculations
	myformula <- as.formula(print(paste(y,x, sep="~")));
	mylm <- lm(myformula, data=data, ...);
	myLevene <- leveneTest(myformula, data=data, center=median);
	myShapiro <- shapiro.test(data[[y]]);
	standard.deviations <- lapply(split(data[[y]], data[[x]]), sd,na.rm=TRUE);
	means <- lapply(split(data[[y]], data[[x]]), mean,na.rm=TRUE);
	n.length <- lapply(split(data[[y]], data[[x]]), length);
	boxplot(myformula, data=data);
	
	#ggplot2 does not work for now :-(
	#myboxplot <- ggplot(data=data, aes_string(x=x, y=y)) + geom_boxplot(aes_string(fill=x));
	#print(myboxplot);
	
	myanova <- as.data.frame(anova(mylm));
	myanova.ss <- myanova$"Sum Sq";
	Eta <- myanova.ss/(myanova.ss + myanova.ss[length(myanova.ss)])
	#anova output
	output <- list();
	colnames(myanova) <- c("DF", "SumSq", "MeanSq", "F-value", "P");
	standard.deviations <- lapply(standard.deviations, unname);
	means <- lapply(means, unname);
	n.length <- lapply(n.length, unname);
	pairwise <- pairwise.t.test(data[[y]], data[[x]], p.adj="bonf");
	output$terms <- apply(myanova,1,as.list);
	output$Eta <- Eta[1]
	output$levene <- unclass(myLevene);
	output$shapiro <- myShapiro$p.value;
	output$pairwise <- unclass(pairwise);
	output$means <- means;
	output$standard.deviations <- standard.deviations;
	output$n <- n.length;
	#Better to do a regression for this:
	#coefficients output
	#output$coef <- data.frame(
	#	coef=coef(mylm), 
	#	assign=mylm$assign, 
	#	terms=c("(Intercept)",attr(terms(mylm),"term.labels")[mylm$assign])
	#);
	
	return(output);
}


 #EXAMPLE: 
#library(car)
#library(foreign)
#SPSS <- read.spss("//Statisticssolut/shared docs/Results/ANOVA/ANOVA data.sav")
#oneway.anova("Species", "Sepal.Width", iris)
#oneway.anova("Race","Test2", SPSS)
#library(opencpu.encode);
#cat(asJSON(oneway.anova("Race","Test2", SPSS)));

#
#

#setwd("C:/Documents and Settings/Administrator/Desktop")
#Data<-read.table("Rori Test2.csv", header = TRUE)
#oneway.anova(Ethnicity2,Math1,Data))
