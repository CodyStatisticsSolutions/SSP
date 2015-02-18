#' Exploratory Factor Analysis
#' Required Libraries : psych, parallel, MASS, GPArotation, MVN, lattice
#' 
#' @Perform  Exploratory Factor Analysis
#' @param df - data frame to derive correlation matrix
#' @param rotate - type of rotation to transform factor loading, defaulted to "varimax"
#' @param user.provide.nfact - Is user providing number of factors to be extracted, defaulted to FALSE i.e. the procedure with find the optimal number
#'                             based on what method is used (see nfact.select.method)
#' @param nfact.from.user - Number of factors provided by user                            
#' @param nfact.select.method - methods to select number of factors, defaulted to "parallel" i.e. using parallel analysis
#' @param fm - methods to extract factors, defaulted to maximum likelihood ("ml")
#' @param use - what observations to be used to calculate correlation matrix, defaulted to using complete observations ("complete.obs") 
#' @param ... - additional arguments to be passed to sem 
#' @return output required parameters and statistics

#options(digits=4)

#library(psych)
#library(parallel)
#library(MASS)
#library(GPArotation)
#library(MVN)
#library(lattice)
#library(jsonlite)

########################################################
# Codes for figuring out the optimal number of factors #
########################################################

# 1) Comparison data (CD) approach by Ruscio and Roche (2012) #
# Psychol Assess. 2012 Jun;24(2):282-92. doi: 10.1037/a0025697. Epub 2011 Oct 3. #
# load function fa_cd(), fa.parallel1() & scree1() # 
#source("E:\\Statistics Solutions\\Functions\\SS005 - Exploratory Factor Analysis\\cd_parallel.R")

# 2) Parallel analysis & scree plot #
# Horn, John (1965) A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179-185.
# Extracting factors until the eigen values of the real data are less than # 
# the corresponding eigen values of a random data set of the same size #  
fa.parallel.mod <- function(...) {
  do.call("fa.parallel1", list(...)) 
}

# 3) Wayne Velicer's Minimum Average Partial (MAP) criterion #
# Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. 
# Psychometrika, 41, 321-327. doi:10.1007/BF02293557 #
# http://flash.lakeheadu.ca/~jljamies/Nfactors.pdf #
VSS.mod <- function(...) {
  do.call("VSS", list(...)) 
}

# Conduct the actual factor analysis #
fa.mod <- function(...) {
  do.call("fa", list(...)) 
}

########################################################################################################################################
##################################################### Main Function efa() ##############################################################
########################################################################################################################################

efa <- function(data, vars, 
                rotate = c("varimax","oblimin","promax","none"), 
                fm = c("ml","minres","pa"),
                user.provide.nfact = FALSE, nfact.from.user = 1, 
                nfact.select.method = c("cd","parallel","map"), 
                use = c("complete.obs","pairwise"), 
                warnings=FALSE, ...) {
  
  df <- data[,vars]
 
  # Check for duplicate columns #
  # 1) Indicate which columns are duplicate and to be removed #
  # 2) Continue running the codes with unique columns #
  dupvar <- which(duplicated(as.matrix(df), MARGIN=2))
   if (length(dupvar)>0) {
    df <- df[,!names(df) %in% names(df)[dupvar]] 
     } 
  
  # Assumption Tests (Univariate & Multivariate outlier(s) will be handled by separate modules) #
  # 1) Univariate and Multivariate Normal #
  # Shapiro-Wilk's test to examine univariate normality - only use non-missing values #
  uvn.test.results <- do.call("rbind.data.frame", lapply(df,shapiro.test))[,c("statistic","p.value")]
  names(uvn.test.results)[1] <- "Shapiro-Wilk Statistic"
  uvn.test.results[,"variable"] <- rownames(uvn.test.results)
  uvn.test.results[,"Normally Distributed?"] <-
    with(uvn.test.results,ifelse(p.value<0.05,"NO","YES"))
  rownames(uvn.test.results) <- NULL
  
  # Mardia's Test to examine multivariate normality #
  mvn.test.results <- mardiaTest(df[complete.cases(df),], cov=TRUE, qqplot=TRUE)
  mvn.test.results <- list(skewness.p.value = mvn.test.results@p.value.skew,
                           kurtosis.p.value = mvn.test.results@p.value.kurt)
  
  # 2) Linearity #
  # Use Scatter Plot Matrices (splom) #
  print(splom(df, main="Scatter Plot Matrix", xlab=""))

  # Determine optimal number of factors - Default to using comparison data (cd) approach #
  # Set the maximum number of factors to be total number minus 1 #
  # Number of factors has to be less than the number of original variables for the FA to be meaningful #
  if (user.provide.nfact==FALSE)
  {  
   if (nfact.select.method[1]=="cd") {
   # 1) Comparison data approach, with plot #
     nfact.select <- fa.cd(df[complete.cases(df),], F.Max=ncol(df)-1, Graph=T)
     } else 
       if (nfact.select.method[1]=="parallel") {
       # 2) Parallel analysis scree plot approach, with plot #
        capture.output(nfact.parallel <- fa.parallel.mod(df, fm=fm[1], fa="fa", use=use[1])) 
         nfact.select <- nfact.parallel$nfact
         } else {
           # 3) Wayne Velicer's Minimum Average Partial (MAP) criterion #
            nfact.select <- which.min(VSS.mod(df, n=ncol(df)-1, rotate=rotate[1], fm=fm[1], plot=FALSE, use=use[1])$map)
      }
  } else nfact.select <- nfact.from.user
  
  # Generate scree plot if nfact.select.method is not "parallel" #
  if (!nfact.select.method[1]=="parallel") scree1(rx=df, use=use[1], pc=FALSE)
  
  # Conduct Exploratory Factor Analysis using the optimal number of factors or provided by user #
  fa.results <- fa.mod(df, max.iter=1e8, nfactors=nfact.select, fm=fm[1], rotate=rotate[1], use=use[1], warnings=warnings, ...)
              
  # Output relevant FA components # 
  # Retrieve factor loadings, communality #
  loadings <- as.data.frame(unclass(fa.results$loadings))
  names(loadings) <- paste("factor",1:ncol(loadings))
  # Add communality to the factor loadings #
  loadings[,"communality"] <- fa.results$communality
  loadings[,"varname"] <- row.names(loadings)
  row.names(loadings) <- NULL
  
  # Fit Statistics #
  RMSEA <- fa.results$RMSEA
  TLI <- fa.results$TLI
  
  # Calculate SS Loadings and related stats #
  square = function(x) x * x
  # Sum of Square Loadings #
  SS.loadings <- apply(as.data.frame(loadings[,1:nfact.select]),2,function(x) sum(square(x)))
  # Proportion of Variance Explained out of total variance #
  Var.proportion <- SS.loadings / nrow(loadings)
  # Cumulative Variance Explained out of total variance #
  Var.cumulative <- cumsum(Var.proportion)
  # Proportion of variance Explained out of all factors (sum to 100%) #
  Var.proportion.factor <- SS.loadings / sum(SS.loadings)
  # Cumulative Variance Explained out of all factors (sum to 100%) # 
  Var.Cumulative.factor <- cumsum(Var.proportion.factor)
  SS.loadings <- as.data.frame(rbind(SS.loadings,Var.proportion,Var.cumulative,Var.proportion.factor,Var.Cumulative.factor))
  names(SS.loadings) <- paste("factor",1:ncol(SS.loadings))
  SS.loadings[,"attributes"] <- row.names(SS.loadings)
  row.names(SS.loadings) <- NULL
  
  # Outputs to json #
  return(list(uvn.test.results = uvn.test.results , 
              mvn.test.results = mvn.test.results, 
              nfact.select = nfact.select, 
              loadings = loadings,
              RMSEA = RMSEA,
              TLI = TLI,
              SS.loadings = SS.loadings
  ))  
    
}



########################################################################################################################################
################################################## End of Main Function efa() ##########################################################
########################################################################################################################################

# Testing #

# Using test data #
#testdata <- read.csv("E:\\Statistics Solutions\\testdata.csv")
#data <- read.csv(file=file.choose())
# Testing without missing data #
#data.nm <- testdata[,c("NormalDependent","Normal","Outliers1","Paired4","IVModerate")]
# Testing with missing data #
#data.m <- data[,c("NormalDependent","Normal","Outliers1","Paired4","Paired2Different","IVModerate","NoCorrelate","Outliers2")]

# run test codes #
#.rotate <- c("varimax","oblimin","promax","none") 
#.fm <- c("ml","minres","pa") 
#.nfact.select.method <- c("cd","parallel","map") 
#.use <- c("complete.obs","pairwise")

# no missing data #
#for (.f in .fm) {
# for (.r in .rotate) {
#  for (.n in .nfact.select.method) {
#   for (.u in .use) {
#     print(c(.r,.n,.u,.f))
#     #efa(data.nm,nfact.select.method=.n,use=.use,rotate=.r,fm=.f)
#     efa(data.nm,nfact.select.method=.n,use=.use,rotate=.r,fm=.f,user.provide.nfact=TRUE,nfact.from.user=1)
#   }
#  } 
# }
#}

# with missing data #
#for (.f in .fm) {
#  for (.r in .rotate) {
#    for (.n in .nfact.select.method) {
#      for (.u in .use) {
#        print(c(.r,.n,.u,.f))
#        #efa(data.m,nfact.select.method=.n,use=.use,rotate=.r,fm=.f,user.provide.nfact = TRUE, nfact.from.user = 1)
#        efa(data, vars=c("NormalDependent","Normal","Outliers1","Paired4","Paired2Different","IVModerate","NoCorrelate","Outliers2"),nfact.select.method=.n,use=.use,rotate=.r,fm=.f)
#      }
#    } 
# }
#}

#a <- efa(data, vars=c("NormalDependent","Normal","Outliers1","Paired4","Paired2Different","IVModerate","NoCorrelate","Outliers2"),nfact.select.method="cd",use="complete.obs",rotate="varimax",fm="ml") # Warnings but work fine #

#toJSON(a)

#cat(a)
#cat(asJSON(a))
#library(opencpu.encode)

#a

#efa(data, vars=c("NormalDependent","Normal","Outliers1","Paired4","Paired2Different","IVModerate","NoCorrelate","Outliers2"),nfact.select.method="map",use="complete.obs",rotate="oblimin",fm="ml") # Warnings but work fine #







