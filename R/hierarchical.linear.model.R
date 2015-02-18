#' Fit hierarchical models
#' Required Libraries : nlme
#' 
#' @Perform  hierarchical modelling
#' @param y - continuous response variable
#' @param fixedx and randomx - list of fixed and random effects of the alternative model to be inputted into the model #
#' @param fixedx.null and randomx.null - list of fixed and random effects of the null model to be inputted into the model #
#'                                     - Variables are put as a list of variables without quotes and separated by comma  #
#'                                     - e.g. fixedx=list(x1,x2,x3) and randomx=list(r1,r2,r3) #
#' @param grpvar - Grouping variable
#' @param data - data frame that contains the data
#' @param method - either Maximum Likelihood ("ML") as default, or Restricted Maximum Likelihood ("REML")
#' @param form - What QQ normal plot to draw, default to standardized residuals
#' @param ... - additional arguments to be passed to lme.mod(), a wrapper for lme()

#' @return output - parameters and statistics

#options(digits=4)

#library(nlme)

# Allow existing functions to take external parameters within the function lme2 #
lme.mod <- function(...) {
  do.call("lme", list(...)) 
}

# REML yields better random effect estimates, but for model comparsion we set ML as the method default #
lme2 <- function(y, 
                 fixedx=list(1), fixedx.null=list(1), 
                 randomx=list(1), randomx.null=list(1),
                 grpvar, data, method="ML", form=~resid(.,type="p"), ...) {
  
  # Combine the fixed effects with response into a formula expression #  
  y <- deparse(substitute(y))    
  response <- paste0(y,"~")
  
  fixedx <- paste(substitute(fixedx))[-1]
  response.fixedx <- paste(fixedx,collapse="+")
  response.fixedx <- as.formula(paste0(response,response.fixedx))
  
  fixedx.null <- paste(substitute(fixedx.null))[-1]
  response.fixedx.null <- paste(fixedx.null,collapse="+") 
  response.fixedx.null <- as.formula(paste0(response,response.fixedx.null))
  
  grpvar <- deparse(substitute(grpvar))
  
  randomx <- paste(substitute(randomx))[-1]
  randomx.grpvar <- paste(randomx,collapse="+")
  randomx.grpvar <- as.formula(paste0("~",randomx.grpvar,"|",grpvar))
  
  randomx.null <- paste(substitute(randomx.null))[-1]
  randomx.null.grpvar <- paste(randomx.null,collapse="+")
  randomx.null.grpvar <- as.formula(paste0("~",randomx.null.grpvar,"|",grpvar)) 
  
  # Select complete cases for modelling #
  allxy <- unique(c(y,fixedx,fixedx.null,randomx,randomx.null,grpvar))
  allxy <- allxy[-which(allxy=="1")]
  data <- data[complete.cases(data[,allxy]),]
  
  # Set up optimizer to increase # of iteration to convergence #
  ctrl <- lmeControl(opt='optim')
  
  # Fitting Alternative Model # 
  results.fitted <- lme.mod(fixed=response.fixedx, random=randomx.grpvar, data=data, control=ctrl, method=method, ...)
  
  # Fitting Null Model #
  results.null <- lme.mod(fixed=response.fixedx.null, random=randomx.null.grpvar, data=data, control=ctrl, method=method, ...)
  
  # Fixed Effect Coefficients #
  fixedx.coeffs <- as.data.frame(summary(results.fitted)$tTable) 
  names(fixedx.coeffs) <- c("B","SE","DF","t","p.value")
  
  # Variance Component of random effects and unexplained within Subject residuals #
  randomx.varcorr <- data.frame(Variance=as.matrix(nlme::VarCorr(results.fitted))[,1],
                                stringsAsFactors=FALSE)
  randomx.varcorr[,"Variance"] <- as.numeric(randomx.varcorr[,"Variance"])
  randomx.varcorr[,"% of Total Variance"] <- round(100*randomx.varcorr[,"Variance"]/sum(randomx.varcorr[,"Variance"]),2)
  
  # Model Fit Statistics #
  AIC.null <- AIC(results.null)
  AIC.fitted <- AIC(results.fitted)
  BIC.null <- BIC(results.null)
  BIC.fitted <- BIC(results.fitted)
  anova.fitted <- suppressWarnings(anova(results.fitted,results.null))
  anova.fitted.LR <- anova.fitted$L.Ratio[2]
  anova.fitted.p.value <- anova.fitted$"p-value"[2]
  model.fit.stat <- c("AIC.null","AIC.fitted","BIC.null","BIC.fitted","anova.fitted.LR","anova.fitted.p.value")
  model.fit.value <- c(AIC.null,AIC.fitted,BIC.null,BIC.fitted,anova.fitted.LR,anova.fitted.p.value)
  
  if (method=="REML") 
  {reminder <- warning("REML anova comparisons are not meaningful if fitted models have different fixed effects", call. = FALSE)}
  else 
  {reminder <- warning("For Likelihood Ratio Test make sure the null model is nested within the alterative",call. = FALSE)}
  
  # Homogeneity of Variance with each group #
  residuals.df <- data.frame(res=residuals(results.fitted),group=names(residuals(results.fitted)))
  
  # Bartlett's test of equal variance #
  bartlett <- with(residuals.df,bartlett.test(res,group))
  bartlett <- paste(bartlett$method,"p-value=",round(bartlett$p.value,4))
  
  # Attach Bartlett's test results to the plot #
  boxplot(res~group, data=residuals.df, ylab="Residuals", xlab=grpvar, 
          main=paste("Residual Boxplot by",grpvar,"\n",bartlett))
  abline(h=0,col="red")
  
  # Investigate standardized residuals #
  # normal plot of the standardized residuals evaluated at the innermost level of nesting
  plot(qqnorm(y=results.fitted, form=form, abline=c(0,1)))
  
  return(list(fixedx.coeffs=fixedx.coeffs,
              randomx.varcorr=randomx.varcorr,
              model.fit.stat=model.fit.stat,
              model.fit.value=model.fit.value,             
              reminder=reminder
  ))
  
}
  
#####################################################################################

# test data #

#data(Orthodont,package="nlme")
#Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
#lme2(y=distance,fixedx=list(age),randomx=list(age,nsex),randomx.null=list(age,nsex),data=Orthodont,grpvar=Subject)

#data3 <- group.center("Correlate", "Group3", "Correlate.Group", data2)

#lme2(y=repdat, fixedx=list(Group3, Correlate), randomx=list(Correlate), data=data2, grpvar=contrasts)

#lme2(y=repdat, data=data2, grpvar=contrasts)
#library(jsonlite)
#library(opencpu.encode)
#cat(asJSON(lme2(y=repdat, fixedx=list(Group3, Correlate.Mean), randomx=list(Correlate.Mean), data=data4, grpvar=contrasts)))
#toJSON(lme2(y=repdat, fixedx=list(Group3, Correlate), randomx=list(Correlate), data=data2, grpvar=contrasts))
#####################################################################################




#data5 <- stack.data(dvs=c("Paired1", "Paired2", "Paired3", "Paired4"), ivs=c("Group3", "Paired2Different", "Paired2NotDifferent"), id = "ID", data=data)



#data5 <- grand.center("Paired2Different", "Paired2Different.Grand.Mean", data5)
#data5 <- group.center("Paired2NotDifferent", "contrast","Paired2NotDifferent.Group.Mean", data5)

#data5$Paired2Different.Grand.Mean[data5$Paired2Different.Grand.Mean=="NA"] <- 2


#lme2(y=repdat, fixedx=list(Paired2Different.Grand.Mean, Paired2NotDifferent.Group.Mean), data=data5, grpvar=contrasts)