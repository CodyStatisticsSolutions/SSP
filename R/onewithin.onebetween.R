onewithin.onebetween <- function(dvs, iv, id, data){
  
  require(nlme)
  
  
  dataframe <- cbind(data[id], data[dvs], data[iv])
  
  means <- as.data.frame(aggregate(dataframe[dvs], data[iv], mean, na.rm=TRUE))
  sds <- as.data.frame(aggregate(dataframe[dvs], data[iv], sd, na.rm=TRUE))
  
  means2 <- colMeans(dataframe[dvs], na.rm=TRUE)
  
  means <- as.data.frame(t(means))
  colnames(means) <- c(levels(data[[iv]]))
  means <- means[2:nrow(means),]
  
  sds <- as.data.frame(t(sds))
  colnames(sds) <- c(levels(data[[iv]]))
  sds <- sds[2:nrow(sds),]
  
 

  
  shapiros <- sapply(data[dvs], shapiro.test);
  

  
  
  myexpr.df <- make.rm(constant=c(id, iv), repeated=c(dvs), data=dataframe);
  
  
  formula <- eval(parse(text=paste("repdat~contrasts +",iv, "+ contrasts:",iv)));;
  random1 <- eval(as.formula(paste("~1|",id)));
  owobanova <- lme(formula, data=myexpr.df, random=random1, na.action=na.omit);
  
  
  anova <- anova.lme(owobanova, type="sequential", adjustSigma = F);
  
  pairwise1 <- pairwise.t.test(myexpr.df$repdat,myexpr.df$contrasts,p.adj="holm", paired=TRUE, pool.sd=FALSE);
  
  pairwise2 <- pairwise.t.test(myexpr.df$repdat, myexpr.df[[iv]], p.adj="holm");
  
  pairwise3 <- pairwise.t.test(myexpr.df$repdat, myexpr.df$contrasts:myexpr.df[[iv]], p.adj="holm");


 
  
  output <- list();
  output$means <- as.list(as.data.frame(means));
  output$means2 <- as.list(means2);
  output$sds <- as.list(sds);
  output$shapiros <- as.list(as.data.frame(shapiros));
  output$anova2 <- unclass(anova);
  output$pairwise1 <- unclass(pairwise1$p.value)
  output$pairwise2 <- unclass(pairwise2$p.value)
  output$pairwise3 <- unclass(pairwise3$p.value)
  
  return(output)
  
}




#library(nlme)
#data <- read.csv(file=file.choose())


#data2 <- read.csv("C:/Users/Statistics Solutions/Desktop/Examples/Permutation Data.csv")

#onewithin.onebetween(dvs = c("Paired1", "Paired2Different", "Paired2NotDifferent"), iv="Group3",id="ID", data=data);
#onewithin.onebetween(dvs = c("Paired1", "Paired2", "Paired3", "Paired4", "Paired5", "Paired6"), iv="DifferentGroup2",id="ID", data=data2);

#library(opencpu.encode)
#cat(asJSON(onewithin.onebetween(dvs = c("Paired1", "Paired2Different", "Paired2NotDifferent"), iv="Group3",id="ID", data=data2)))

#onewithin.onebetween(dvs=c("Reading1", "Reading2", "Reading3", "Reading4"), iv="Group", id="ID", data=data)

