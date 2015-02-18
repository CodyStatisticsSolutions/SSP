#' Fit Cumulative Link Models
#' Required Libraries : Hmisc and ordinal 
#' 
#' @Perform  Proportional Odds Logistic Regression (Default)
#' @param x - predictors
#' @param y - response, must be an ordered factor defined by user
#' @param data  - data frame
#' @param ... - additional arguments to be passed to polr function which performs the ordinal logistic regression
#' @return output - parameters and statistics

#library(ordinal)
#library(Hmisc)
#library(reshape)

# Regression Modeling Strategies: With Applications to Linear Models, Logistic Regression, and Survival Analysis (Harrell 2001) #
# Main Ordinal Regression function using proportional odds model #

ordinal.reg <- function(formula, ...) {
  
  require(Hmisc)
  require(ordinal)  
  
  # Fitting proportional odds model #
  clm.mod <- function(formula, ...) {
    do.call("clm", list(formula, ...)) 
  }
  
  results.fitted <- clm.mod(formula=formula, ...)
  
  # Null Model #
  
  y <- all.vars(formula)[1]
  formula.null <- as.formula(paste(y," ~ 1"))
  results.null <- clm.mod(formula=formula.null, ...) 
  
  # Proportional Odds assumption plot #
    ordinal.assumption <- function(formula, data, weights=NULL) {
    
    # For weighted data #
    if(!is.null(weights)) 
    { 
      wtvar <- deparse(substitute(weights))
      data <- untable(data[,!names(data) %in% wtvar],data[,wtvar])
    }
    
    resp <- all.vars(formula)[1]
    resp.levels <- levels(data[,resp])
    resp.int <- as.integer(factor(resp.levels,levels=resp.levels))
    check.logit <- paste("check.logit <- function(resp) ",
                         paste0("c(",paste("'resp<=",resp.levels,"'","=qlogis(mean(resp<=",  resp.int,"))",sep="",collapse=","),")"))
    eval(parse(text=check.logit))
    
    # Compare individual logits #
    resp.num <- as.integer(data[,resp])
    formula.check <- as.formula(paste("resp.num ~",paste(all.vars(formula)[-1],collapse="+")))
    s <- summary(formula.check,data=data, fun=check.logit)
    plot(s, which = 1:length(resp.levels), pch = 1:length(resp.levels), 
         xlab = "logit", vnames = "names", xlim = range(c(s[,-1])[which(s[,-1]<Inf & s[,-1]>-Inf)]), 
         main = "Test of Proportional Odds Assumption", width.factor = 2)
    legend("topright", 
           legend = paste0(resp,"<=",resp.levels[1:(length(resp.levels)-1)]), 
           pch = 1:(length(resp.levels)-1),
           cex=0.65,
           ncol=(length(resp.levels)-1))
  }
  
  
  ordinal.assumption(formula=formula,...) 
  
  # Model fit statistics by comparing to null model #
  
  deviance <- 2*(logLik(results.fitted)-logLik(results.null))[1]
  degrees.of.freedom <- df.residual(results.null) - df.residual(results.fitted)
  goodness.of.fit.p <- 1-pchisq(deviance,degrees.of.freedom)
  AIC.null <- AIC(results.null)
  AIC.fitted <- AIC(results.fitted) 
  BIC.null <- BIC(results.null)
  BIC.fitted <- BIC(results.fitted)
  model.fit <- c(deviance=deviance,
                 degrees.of.freedom=as.integer(degrees.of.freedom), 
                 goodness.of.fit.p=goodness.of.fit.p,
                 AIC.null=AIC.null,AIC.fitted=AIC.fitted, 
                 BIC.null=BIC.null,BIC.fitted=BIC.fitted)
  
  # Parameter Estimates and OR including 95% CI (ORLL, ORUL) #
  
  param <- as.data.frame(summary(results.fitted)$coefficients)
  names(param) <- c("B","SE","Z","p.value")
  param[,"OR"] <- exp(param[,"B"])
  param[,"ORLL"] <- exp(param[,"B"]-1.96*param[,"SE"])
  param[,"ORUL"] <- exp(param[,"B"]+1.96*param[,"SE"]) 
  
  # Output relevant statistics #
  
  output <- list(model.fit=model.fit, param=param)
  return(output)
  
}

#data <- read.csv(file=file.choose())

#data$TypicalLikert3 <- ordered(data$TypicalLikert3)

#ordinal.reg(TypicalLikert3~NotDifferentGroup3 + ChiGroup2Large,data=data)


#ordinal.reg(Sat~Infl,data=housing,weights=Freq)
#library(opencpu.encode)
#cat(asJSON(ordinal.reg(TypicalLikert3~Paired1 + Paired2Different,data=data)))
