#' Fit Structural Equation Modelling
#' Required Libraries : sem
#' 
#' @Perform  Structural Equation Modelling (Default)
#' @param model - structural equations
#' @param data - data frame to derive correlation matrix
#' @param ... - additional arguments to be passed to sem 
#' @return output - parameters and statistics

#options(digits=4)

#library(sem)
#library(MVN) # for multivariate normal assumption using Mardia test #
#library(stringr)
#library(semPlot) # Plot path diagram for SEM models #

# Allow clm to take external parameters within the function ordinal.reg #
# non-normed fit index (NNFI) is the same as Tucker-Lewis index (TLI) #



makeSemModel <- function(modelSpec, semdata){
  tmp <- tempfile()
  writeLines(modelSpec, tmp)
  semmodel <- specifyEquations(file=tmp)
  sem2(semmodel, data=semdata, what="est", layout="tree2")
}


# Modified Modification Index function #

modIndices.modified <- function(x, n.largest.a=5, n.largest.p=5, ...){
  printa <- printp <- A.path <- P.path <- NULL
  mod.A <- as.vector(x$mod.A)
  .names <- rownames(x$mod.A)
  names(mod.A) <- outer(.names, .names, paste, sep="<-")
  if (n.largest.a>0) {
    printa <- (rev(sort(mod.A))[1:n.largest.a])
    A.path <- names(printa) 
  }
  mod.P <- as.vector(x$mod.P)
  names(mod.P) <- outer(.names, .names, paste, sep="<->")
  if (n.largest.p>0) {
    printp <- (rev(sort(mod.P[lower.tri(x$mod.P, diag=TRUE)]))[1:n.largest.p])
    P.path <- names(printp)  
  }
  output <- list(A.path=A.path,mod.indices.A.Matrix=printa,P.path=P.path,mod.indices.P.Matrix=printp)
  return(output)
}


sem.mod <- function(...) {
  do.call("sem", list(...)) 
}


semPaths.mod <- function(...) {
  do.call("semPaths", list(...)) 
}

sem2 <- function(model, data, ...) {
  
  require(sem)
  require(semPlot)
  
  
  # Select the variables actually used in SEM #
  findvar <- str_locate_all(model[,1]," ")
  lhs <- unlist(lapply(findvar,function(x) x[1,1])) # variables on the LHS #
  lhsvar <- substr(model[,1],1,lhs-1)
  rhs <- unlist(lapply(findvar,function(x) x[2,1])) # variables on the RHS #
  rhsvar <- substr(model[,1],rhs+1,str_length(model[,1]))
  allvar <- sort(unique(c(lhsvar,rhsvar)))
  data <- data[,intersect(names(data),allvar)]  
  
  # Actual Sem begins #
  sem.results <- sem.mod(model=model, data=data, ...)
  sem.summary <- summary(sem.results, digits=getOption("digits"), 
                         conf.level=.90, robust=FALSE,  
                         fit.indices=c("RMSEA","NNFI", "CFI"))
  
  # Fit Statistics #
  .chisq <- round(sem.summary$adj.obj$chisq,4)
  .df <- sem.summary$adj.obj$df
  .p <- round(sem.summary$adj.obj$p.old,4)
  .CFI <- round(sem.summary$CFI,4)
  .TLI <- round(sem.summary$NNFI,4)
  .RMSEA <- round(sem.summary$RMSEA[1],4)
  .RMSEA.90CI <- paste0("(",paste(round(sem.summary$RMSEA[2:3],4),collapse=","),")")
  
  fit.stat <- list(chi.sq=.chisq,
                   df=.df,
                   p.value=.p,
                   CFI=.CFI,
                   TLI=.TLI,
                   RMSEA=.RMSEA,
                   RMSEA.90CI=.RMSEA.90CI)
  
  
  # Modification Indices >=10 will be displayed #
  # The 'A' matrix represents regression paths:  x:y represents y->x .  
  # The 'P' matrix represents covariance paths: x:y represents x<->y .
  
  modind <- modIndices(sem.results)
  A.gt.10 <- length(which(modind[[1]]>=10))
  P.gt.10 <- length(which(modind[[2]]>=10)) 
  modind.greater.than.10 <- modIndices.modified(modind, n.largest.a=A.gt.10, n.largest.p=P.gt.10)
  
  # Unstandardized Parameter Estimates #
  unstd.param.est <- data.frame(sem.summary$coeff,stringsAsFactors=FALSE)
  names(unstd.param.est) <- c("Unstd. Est.","StdError","z","p","path")
  unstd.param.est[,"param"] <- rownames(unstd.param.est)
  rownames(unstd.param.est) <- NULL
  
  # Standardized Parameter Estimates #
  std.param.est <- data.frame(standardizedCoefficients(sem.results),
                              stringsAsFactors=FALSE)
  names(std.param.est) <- c("param","Std. Est.","path")
  
  # Combine into one parameter estimate data frame #
  param.est <- merge(unstd.param.est,std.param.est,by=c("param","path"),all=TRUE)
  param.est <- param.est[,c("param","path","Unstd. Est.","Std. Est.","StdError","z","p")]
  
  # Shapiro-Wilk's test to examine univariate normality #
  uvn.test.results <- as.data.frame(t(sapply(data,shapiro.test)))[,c("statistic","p.value")]
  names(uvn.test.results)[1] <- "Shapiro-Wilk Statistic"
  uvn.test.results[,"variable"] <- rownames(uvn.test.results)
  uvn.test.results[,"Normally Distributed?"] <-
    with(uvn.test.results,ifelse(p.value<0.05,"NO","YES"))
  rownames(uvn.test.results) <- NULL
  
  # Mardia's Test to examine multivariate normality #
  mvn.test.results <- mardiaTest(data, cov = TRUE, qqplot = TRUE)
  
  # Graphing SEM paths #
  sem.plots <- semPaths.mod(
    sem.results,
    residuals=FALSE,
    intercepts=FALSE,
    nCharNodes=0,
    nCharEdges=0,
    asize=2,sizeMan=4,
    label.prop=1.2,
    posCol="black",
    DoNotPlot=FALSE,...) 
  
  # Output results #
  output <- list(fit.stat=fit.stat, 
                 param.est=param.est,
                 modind.greater.than.10=modind.greater.than.10,
                 uvn.test.results=unlist(uvn.test.results),
                 mvn.test.results=mvn.test.results,
                 sem.plots=sem.plots
  )
  return(output)
  
}

#####################################################################################

# test data #

#modelbollen <- "
#y1 = 1*Demo60
#y2 = l1.2*Demo60
#y3 = l1.3*Demo60
#y4 = l1.4*Demo60
#y5 = 1*Demo65
#y6 = l2.2*Demo65
#y7 = l2.3*Demo65
#y8 = l2.4*Demo65
#x1 = 1*Indust
#x2 = l3.2*Indust
#x3 = l3.3*Indust
#c(y1,y5)=c1
#c(y2,y4)=c2
#c(y2,y6)=c3
#c(y3,y7)=c4
#c(y4,y8)=c5
#c(y6,y8)=c6
#Demo60 = b1.1*Indust
#Demo65 = b2.1*Indust + b2.2*Demo60
#V(Indust) = v1
#V(y1)=v2
#V(y2)=v3
#V(y3)=v4
#V(y4)=v5
#V(y5)=v6
#V(y6)=v7
#V(y7)=v8
#V(y8)=v9
#V(x1)=v10
#V(x2)=V11
#V(x3)=v12
#V(Demo60) = v13
#V(Demo65) = v14
#"

#makeSemModel(modelbollen, semdata=Bollen)



#sem.bollen
#library(opencpu.encode)
#cat(asJSON(sem.bollen))

#####################################################################################





#semdata <- read.csv(file=file.choose())


#sem.model <- "
#CWD01_I = 1*CWD
#CWD02_I = l1.2*CWD
#CWD03_I = l1.3*CWD
#CWD04_I = l1.4*CWD
#CWD05_I = l1.5*CWD
#AWDI11_I = 1*AWD
#AWDI12_I = l2.2*AWD
#AWDI13_I = l2.3*AWD
#AWDI14_I = l2.4*AWD
#AWDI15_I = l2.5*AWD
#RDD21_I = 1*RDD
#RDD22_I = l3.2*RDD
#RDD23_I = l3.3*RDD
#RDD24_I = l3.4*RDD
#RDD25_I = l3.5*RDD
#Oscillation = b1.1*CWD + b1.2*AWD + b1.3*RDD
#ComplicatedGrief = b2.1*Oscillation + b2.2*CWD + b2.3*AWD + b2.4*RDD
#CWD = b3.1*AWD + b3.2*RDD
#V(CWD01_I) =v1
#V(CWD02_I) =v2
#V(CWD03_I) =v3
#V(CWD04_I) =v4
#V(CWD05_I) =v5
#V(AWDI11_I) = v6
#V(AWDI12_I) = v7
#V(AWDI13_I) = v8
#V(AWDI14_I) = v9
#V(AWDI15_I) = v10
#V(RDD21_I) = v11
#V(RDD22_I) = v12
#V(RDD23_I) = v13
#V(RDD24_I) = v14
#V(RDD25_I) = v15
#V(Oscillation) = v16
#V(ComplicatedGrief) = v17
#V(CWD)=v18
#V(AWD)=v19
#V(RDD)=v20
#"

#makeSemModel(model=sem.model, semdata)



