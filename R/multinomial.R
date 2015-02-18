#' Fit Multinomial Logistic Regression
#' Required Libraries : nnet, plyr, reshape2 
#' 
#' @Perform  Multinomial Odds Logistic Regression 
#' @param formula - a formula expression as for regression models, of the form response ~ predictors
#' @param data  - data frame
#' @param ref - reference group 
#' @param ... - additional arguments to be passed to function mlr 
#' @return output - parameters and statistics

#options(digits=4)

#library(nnet)
#library(plyr)
#library(reshape2)
#library(jsonlite)

# Allow multinom to take external parameters within the function mlr #

multinom.mod <- function(...) {
  do.call("multinom", list(...)) 
}

# Main Multinomial Logistic Regression codes #
# R uses lowest factor level as default reference #
# Use relevel if user wants a difference reference group #

mlr <- function(formula, data, weights=NULL, ref="", trace=FALSE, ...) {
  
  # Select complete cases for modelling #
  data <- data[complete.cases(data[,all.vars(formula)[-1]]),]
  
  # Response variable has to be a factor - transform response matrix into long format prior to modelling # 
  # Relevel the response variable if required #
  y <- all.vars(formula)[1]
  if (ref != "")
  {
   data[,y] <- relevel(data[,y],ref=ref)
  }
  
  # Determine weights - default to 1 #
  wts <- deparse(substitute(weights))
  if (wts == "NULL")
   {wts <- rep(1,nrow(data))}
  else
   {wts <- data[,wts]}
  
  # Fitting proportional odds model #  
  results.fitted <- multinom.mod(formula=formula, data=data, trace=trace, weights=wts, ...)
  
  # Null Model #
  formula.null <- as.formula(paste(y," ~ 1"))
  results.null <- multinom.mod(formula=formula.null, data=data, trace=trace, weights=wts, ...) 
  
  # Model fit statistics by comparing to null model #
  dev <- deviance(results.null) - deviance(results.fitted)
  degrees.of.freedom <- results.fitted$edf - results.null$edf
  goodness.of.fit.p <- 1-pchisq(dev,degrees.of.freedom)
  AIC.null <- AIC(results.null)
  AIC.fitted <- AIC(results.fitted) 
  BIC.null <- BIC(results.null)
  BIC.fitted <- BIC(results.fitted)
  
  model.fit.stat <- c("deviance","degrees.of.freedom","goodness.of.fit.p","AIC.null","AIC.fitted","BIC.null","BIC.fitted")
  model.fit.value <- c(dev,as.integer(degrees.of.freedom),goodness.of.fit.p,AIC.null,AIC.fitted,BIC.null,BIC.fitted)
  
  # Parameter Estimates and Relative Risk (RR) including 95% CI (RRLL, RRUL) #
  param <- as.data.frame(summary(results.fitted)$coefficients)
  param[,"Response.Category"] <- row.names(param)
  param <- melt(param,id="Response.Category")
  names(param)[2:3] <- c("Parameter","B")
  
  stderr <- as.data.frame(summary(results.fitted)$standard.errors)
  stderr[,"Response.Category"] <- row.names(stderr)
  stderr <- melt(stderr,id="Response.Category")
  
  param[,"SE"] <- stderr[,"value"]
  param[,"Z"] <- with(param,B/SE)
  param[,"p.value"] <- with(param,2*(1-pnorm(abs(Z)))) # 2-tailed z test #
  
  # The coefficients in the multinom call are log odds. the exp(coef()) call gets relative risk. #
  param[,"RR"] <- exp(param[,"B"])
  param[,"RRLL"] <- exp(param[,"B"]-1.96*param[,"SE"])
  param[,"RRUL"] <- exp(param[,"B"]+1.96*param[,"SE"]) 
  param <- arrange(param,Response.Category)
  row.names(param) <- NULL
    
  # Output relevant statistics #
  output <- list(model.fit.stat=model.fit.stat, model.fit.value=model.fit.value, param=param)
  return(output)
  
}


#################################################################################################

# Test #

#library(foreign)
#ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
#ml$prog2 <- relevel(ml$prog, ref = "academic")
#test <- multinom(prog2 ~ ses + write, data = ml)
#test0 <- multinom(prog2 ~ 1, data = ml)

#a <- mlr(prog ~ ses + write, data = ml, ref=academic)
#toJSON(a)

# Sample data test #
#library(Hmisc)
#testdata <- read.csv("E:\\Statistics Solutions\\Functions\\SS004 - Multinomial Logistic Regression\\testdata.csv")

#testdata[,"TypicalLikert1"] <- factor(testdata[,"TypicalLikert1"])
#testdata[,"TypicalLikert2"] <- factor(testdata[,"TypicalLikert2"])
#testdata[,"BigGroup"] <- factor(testdata[,"BigGroup"])

# Unweighted #
#a <- mlr(TypicalLikert1 ~ NotNormalDependent + Paired2Different + Paired2NotDifferent +
#                          Outliers1 + Outliers2 + TypicalLikert2 + Correlate +  
#                          PositiveCorrelate + ChiGroup2Small , 
#         data=testdata)
#toJSON(a)

# Weighted #
#b <- mlr(TypicalLikert1 ~ NotNormalDependent + Paired2Different + Paired2NotDifferent +
#           Outliers1 + Outliers2 + TypicalLikert2 + Correlate +  
#           PositiveCorrelate + ChiGroup2Small , 
#         data=testdata,
#         weights = NormalDependent,
#         ref="4")
#toJSON(b)





