kruskal.wallis <- function (x, y, data, ...){
  
  
  formula <- as.formula(paste(y, "~", x));
  test <- kruskal.test(formula, data=data);
  
  
  r <- rank(data[[y]]);
  data2 <- data.frame(r, data[[x]]);
  ranks <- lapply(split(r, data[[x]]), mean, na.rm=TRUE);
  
  pairwise <- kruskalmc(formula, data)
  pairwise2 <- cbind(row.names(pairwise$dif.com), pairwise$dif.com[,3])
  pairwise2 <- as.data.frame(pairwise2)
  
  
  
  boxplot(formula, data = data, main=paste("Scores on", y, "by", x), xlab = x, ylab = y);
  
  pairwise2$V1 <- paste(pairwise2$V1, sep="")
  pairwise2$V2 <- as.logical(pairwise2$V2)+0
  
  output <- as.list(unclass(test));
  output$ranks <- lapply(ranks, unname);
  output$pairwise <- as.list(pairwise2)
return(output);
  
  
}

#data <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")
#kruskal.wallis("Group3", "NotNormalDependent", data)
#kruskal.wallis("ChiSquareNotSig", "NotNormalDependent", data)

#library(pgirmess)

#library(opencpu.encode)
#cat(asJSON(kruskal.wallis("Group3", "NotNormalDependent", data)))

