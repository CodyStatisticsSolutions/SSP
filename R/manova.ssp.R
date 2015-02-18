#library(foreign)
#library(car)
#library(heplots)
#library(psych)


manova.ssp <- function(dvs, x, data, ...){
  #remove NA
  data <- cbind(data[dvs], data[x])
  data <- na.omit(data)
  matrix <- as.matrix(data[dvs])
  
  numdv <- ncol(matrix)
  numgr <- nlevels(data[[x]])
  grlvls <- levels(data[[x]])
  
  
  
  
  #calculations
  myformula <- as.formula(paste("matrix~" ,x));
  mylm <- manova(myformula, data=data);
  
  #myLevene <- leveneTest(myformula, data=data, center=mean);
  #myShapiro <- lapply(data[dvs], shapiro.test);
  
  levenes <- list()
  for (i in 1:numdv){
    levenes[i] <- list(leveneTest(data[,i]~data[[x]])$`Pr(>F)`[1])
  }
  
  shapiros <- list()
  for (i in 1:numdv){
    shapiros[i] <- list(shapiro.test(data[,i])$p.value)
  }
  
  
  
  descriptives <- describeBy(data[dvs], data[x], mat=TRUE);
  n.length <- lapply(split(data[dvs], data[x]), nrow);
  n.length <- sapply(n.length, unname);
  
  descriptives <- cbind(descriptives$group1, descriptives$var, descriptives$mean, descriptives$sd)
  descriptives <- as.data.frame(descriptives)
  descriptives$V2 <- factor(descriptives$V2, levels=c(1:numdv), labels=dvs)
  descriptives$V1 <- factor(descriptives$V1, levels=c(1:numgr), labels=grlvls)
  colnames(descriptives) <- c("IV", "DVs", "Mean", "SD")


  myanova <- Anova(mylm, anova=TRUE)
  anovas <- summary.aov(mylm)
  eta <- etasq(Anova(mylm), anova=TRUE)



  
  ncolpair <- (numdv*(numdv+1))/2  
  pairwise <- list()
  for (i in 1:numdv){
    pairwise[i] <- list(pairwise.t.test(data[,i], data[[x]], p.method="bonferroni")$p.value)
  }  
  names(pairwise) <- dvs

  
  
  #output
  output <- list();  
  MANOVA.F <- anova(mylm)$`approx F`[2]
  MANOVA.df1 <- anova(mylm)$`num Df`[2]
  MANOVA.df2 <- anova(mylm)$`den Df`[2]
  MANOVA.p <- anova(mylm)$`Pr(>F)`[2]
  
  output$shapiro <- shapiros;  
  output$levene <- levenes;
  output$MANOVA.F <- MANOVA.F
  output$MANOVA.df1 <- MANOVA.df1
  output$MANOVA.df2 <- MANOVA.df2
  output$MANOVA.p <- MANOVA.p
  output$ANOVAs <- unclass(anovas)
  output$Eta <- eta$`eta^2`
  output$pairwise <- pairwise;
  output$descriptives <- descriptives;
  output$n <- n.length;
  
  return(output);
}


#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")

#manova.ssp(dvs = c('Paired1', 'Paired2', 'Paired3', 'Paired4', 'Paired5'), x="Group3", data=data);
#manova.ssp(dvs = c('NormalDependent', 'NotNormalDependent', 'Normal'), x="DifferentGroup3", data=data);


#library(opencpu.encode);
#cat(asJSON(manova.ssp(dvs = c('Paired1', 'Paired2Different', 'Paired2NotDifferent'), "DifferentGroup3", data=data)))
#cat(asJSON(manova.ssp(dvs = c('NormalDependent', 'NotNormalDependent', 'Normal'), x="DifferentGroup3", data=data)))
